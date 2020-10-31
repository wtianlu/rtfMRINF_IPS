% Main script for WP3 Neurofeedback
% Online experiment script: reads *.rtp files from Turbo Brainvoyager and
% presents differential PSC with Psychtoolbox
% Lulu 2019.11.19
% ------------------------------------------------
%   IN
% Enter participant # and session, check paths to functions and data 
% directory in section 0. Initialisation
%   OUT
% Saves MATLAB diary and datastructure.mat in participant session folder
% Saves NIFTI, RTP, and DRIN files, requests for eyelink files and 
% transfers to participant folder

%% ----------- %% 0. Initialisation %% ----------- %% 
clc; clearvars; close all;

idP = '1'; % input('Enter participant number [xx]: ','s');
idS = '1'; % input('Enter session number: ','s');
newmaxlvl = 1; % str2double(input('Enter maximum PSC level: ','s'));

% --- Set presentation type
mr8 = ispc; % WinOS at MR8, MacOs in PSI 
bool_disponly = false; % testing phase, set true in simulation mode -> skip TBV, random PSCs
bool_fullscreen = true; % set false for testing with smaller screen

% --- Set code directories
if mr8, cd('C:\Users\fmri\Documents\Lulu\WP3\Code'); 
else, cd('/Users/Lulu/Box Sync/fMRI studies/Neurofeedback/Real-time code');end
addpath([pwd filesep 'functions'])

% --- Initialise participant data structure and folders
ds = initialise_nf(mr8,0,struct(),idP,idS);
diary([ds.dirs.ps filesep ds.idstr '_diary'])
ds.exp.mr8 = mr8; ds.exp.bool_disponly = bool_disponly; ds.exp.bool_fullscreen = bool_fullscreen;

% Get groups, 1 is L-R and 0 is R-L
load part_groups8; ds.group.db.idx = part_groups(str2double(idP)); clear part_groups;
fprintf('%s\tFirst initialisation %s finished\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),ds.idstr)

%% ----------- %% 1. Start Eyelink and record calibration %% ----------- %% 

eye = NF_eyelink('init');
if eye
    disp('Press enter to start Eyelink 9-grid recording');commandwindow; pause;
    
    % -- psychtoolbox screen 
    pt = load_PTBsettings(bool_fullscreen,1);
    Screen('Preference', 'SkipSyncTests', 1);
    [pt.win.w, pt.win.windowRect] = PsychImaging('OpenWindow', mr8*pt.win.screenNumber, pt.stim.colors.dgrey);
    pt = load_PTBsettings(bool_fullscreen,2,pt);
    
    % -- Eyelink recording
    edfFile = ['P' idP 'S' idS 'CA']; el.edfFile = edfFile;
    fprintf('%s\tStart recording EDF %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),edfFile); 
    NF_eyelink('startr',edfFile,pt);    NF_eyelink('calib',edfFile,pt);    NF_eyelink('stopr',edfFile); 
    Screen('closeall');
    fprintf('----------------------\nFinished 9-grid recording\n----------------------\n')
else
    disp(['Eyelink not connected, skip 9-grid calibration (P' idP 'S' idS 'CA)'])
end
fprintf('\n\n')

%% Set up eye tracking for resting-state (run only for session 1 and 3)
edfFile = ['P' idP 'S' idS 'RS2']; el.edfFile = edfFile;
fprintf('%s\tStart recording EDF %s\nWaiting for trigger...\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),edfFile); 
triggertime = GetSecs;
while true
    [keyIsDown, keyTime, keyCode] = KbCheck(pt.input.device);
    if keyCode(pt.input.triggerkey), triggertime = keyTime; break; end
    if keyCode(pt.input.abortkey), disp('>>> Experiment Aborted <<<');return; end
end
fprintf('%s\tTrigger arrived for resting state: %0.3f\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),triggertime); 
NF_eyelink('startr',edfFile,pt);   
%% Stop recording for resting-state
NF_eyelink('stopr',edfFile); 

%% ----------- %% 2. After Struct and Func_test %% ----------- %% 
if ~bool_disponly
    if idS == '1', disp('Start up eyelink recording manually for the resting state scan');end
    input('Press enter if the anatomical/functional test scan has finished')
    fprintf('%s\tStart preparing ROIs\n',datestr(now,'yyyy-mm-dd HH:MM:SS'))
    ds = initialise_nf(mr8,1,ds);

    % -- Transform VOIs from template to native space -> depends on the session
    ds = spm_mni2native(ds);
    fprintf('Check overlays of brwvoi_* on func4D in MRIcron\nPress any key to proceed\n'); pause;

    % -- Transform VOIs from NIFTI to ROI file format
    disp('Transform native ROIs to .roi format');
    nii2tbvroi(ds)
else
    disp('Simulation mode, skipping neuroimaging data and ROIs')
end

%% ----------- %% 3. Prepare file to start  %% ----------- %% 
fprintf('\n\n')
% -- Set current run and prepare experiment variables
run = str2double(input('Enter run nr: ','s'));

% -- check that run hasn't been done before
if exist([ds.dirs.ps filesep 'r' num2str(run)],'dir')~=0
    disp('Run has already been done before! Stop and enter new run. Existing runs:'); dir([ds.dirs.ps filesep 'r*']); pause;
else
    fprintf('Prepare run %d\n',run);
end

% -- initialise run variables
if ~bool_disponly, ds = initialise_nf(mr8,2,ds,'','',run); end % select watch folder and write TBV settings file
if run < 6, istraining = ismember(run,[2 3 4]); else, istraining = str2double(input('test [0] or training [1]: ','s'));end
ds.r(run).istraining = istraining; dirs = ds.dirs; rtp = ds.r(run).rtp; current_files = {}; cur_vol = 0; 

% Option to change the max level at the start of a session (depending on the previous session)
if run == 1 && exist('newmaxlvl','var'), pt.stim.maxlvl = newmaxlvl; end
fprintf('Max feedback level: %.1f\n',pt.stim.maxlvl)

% -- Check that RTP folder is empty
if ~isempty(dir([ds.r(run).dir filesep 'tbv_feedback' filesep '*.rtp']))
    disp(['RTP files present in ' ds.r(run).dir(length(dirs.main)+1:end) ' -> clear contents?'])
    resp = input('[y]/n: ','s');
    if resp == 'n'
        disp('WARNING - Did not delete contents!!!')
    else
        oldfiles = dir([ds.r(run).dir filesep 'tbv_feedback' filesep '*.rtp']);
        oldfiles = fullfile([ds.r(run).dir filesep 'tbv_feedback'],{oldfiles.name})';
        delete(oldfiles{:}); disp('Finished deleting old files')
    end
end

% Check that Eyelink is still connected
res = Eyelink('IsConnected');
if res~=1
    disp('Eyelink is not connected! Re-initialise eyelink...')
    eye = NF_eyelink('init');
    if eye
        disp('Eyelink re-connected!')
    else
        disp('!!! Eyelink connection failed... proceed without automatic eyelink recording!!!')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------- %% 4. Start display  %% ----------- %% 

fprintf('%s\tPress key to start PTB display for run %d\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),run); pause;

pt = load_PTBsettings(bool_fullscreen,1);
Screen('Preference', 'SkipSyncTests', 1);
[pt.win.w, pt.win.windowRect] = PsychImaging('OpenWindow', mr8*pt.win.screenNumber, pt.stim.colors.dgrey);
pt = load_PTBsettings(bool_fullscreen,2,pt);
Screen('TextSize', pt.win.w, pt.win.screenXpixels/24);
pt.win.WindowSize=Screen('WindowSize', pt.win.w);pt.win.Resolution = Screen('Resolution', pt.win.w);

% disp('Set Priority');Priority(1) % set priority to high (WinOS 0=normal, 1=high, 2=real-time)
commandwindow; % go to command window

% --- Show READY screen 
DrawFormattedText(pt.win.w, 'Ready...','center', 'center', pt.stim.colors.white);
[VTs, SOT, FTs] = Screen('Flip', pt.win.w);
ds.r(run).screenflips(1).imgArray = Screen('GetImage',pt.win.w); 
ds.r(run).screenflips(1).timings = [GetSecs, VTs, SOT, FTs];
ds.r(run).screenflips(1).type = 'Ready';

% -- Wait for scan start
fprintf('\n----------------------\n%s\tPTB display presented\nWaiting for trigger...\n',datestr(now,'yyyy-mm-dd HH:MM:SS'))
triggertime = GetSecs;
while true
    [keyIsDown, keyTime, keyCode] = KbCheck(pt.input.device);
    if keyCode(pt.input.triggerkey), triggertime = keyTime; break; end
    if keyCode(pt.input.abortkey), Screen('closeall'); disp('>>> Experiment Aborted <<<');return; end
end
fprintf('%s\tTrigger arrived for run %d: %0.3f\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),run,triggertime); 
ds.r(run).timings.mrtrigger = triggertime;

% -- Check eyelink connection and start recording
if eye
    edfFile = ['P' idP 's' idS 'nf' num2str(run)]; el.edfFile = edfFile;
    try
        fprintf('%s\tStart recording EDF %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),edfFile);
        NF_eyelink('startr',edfFile,pt);
    catch
        fprintf('Eyelink connection error...!\nUse manual recording for %s\n',edfFile);
    end
end

% --- Present fixation for first 10 volumes
Screen(pt.stim.fixation{1},pt.win.w,pt.stim.fixation{2:end}); 
[VTs, SOT, FTs] = Screen('Flip', pt.win.w); 
ds.r(run).screenflips(end+1).imgArray = Screen('GetImage',pt.win.w); %{end+1}
ds.r(run).screenflips(end).timings = [GetSecs, VTs, SOT, FTs];
ds.r(run).screenflips(end).type = 'Fixation';
disp('Cross presented; TBV skipping first ten volumes...')

fprintf('\n--- Checking arrival of DRIN/RTP files ---\n')

% --- Start loop ----------------------
checktime = GetSecs;
while cur_vol<ds.exp.nvol-10
    % Exit loop with escape
    [keyIsDown,keyTime, keyCode, ~] = KbCheck(pt.input.device);
    if keyCode(pt.input.abortkey)
        sca; disp('>>> Experiment Aborted <<<');
        ds.r(run).comments = input('Comments about this run: ','s');
        save([ds.r(run).dir filesep 'ds_abort.mat'],'ds','-v7.3');
        if eye, NF_eyelink('stopr',edfFile); ds.r(run).el = el; end
        return; 
    end
    
    bool_updateScreen = false;
    
    % Check for new files in both DRIN and RTP directories
    if bool_disponly, new_files={};
    else
        dir_contents = [dir([dirs.exportdrin filesep '*.img']); dir([ds.r(run).dir filesep 'tbv_feedback' filesep '*.rtp'])];
        new_files = setdiff({dir_contents.name},current_files);
    end

% -- New files have arrived - save file arrival time
    if ~isempty(new_files)
        fprintf('%.3fs new %s\n',GetSecs-triggertime, new_files{1}); 
    % -- check if new file is RTP file or DRIN
        current_files = {dir_contents.name};
        tmp = strsplit(new_files{1},{'-','.'});
        if isequal(tmp{end},'rtp')
            ds.r(run).timings.rtp = [ds.r(run).timings.rtp;GetSecs];
            cur_vol = str2double(tmp{2}); 
        % -- extract bold signals
            [fileID,errmsg] = fopen([ds.r(run).dir filesep 'tbv_feedback' filesep new_files{1}],'r');
            if isempty(errmsg), [currRTP_tmp,lenRTP] = fscanf(fileID,'%f'); fclose(fileID); else, disp(errmsg);sca;end
            if lenRTP<3,fprintf('!!WARNING!! ROI not loaded in TBV\n'); 
            else, rtp.raw(cur_vol,:) = currRTP_tmp(2:3)'; end
            bool_updateScreen = true;
        else
            ds.r(run).timings.drin = [ds.r(run).timings.drin;GetSecs];
        end
        checktime = GetSecs;
        
% -- New files not arriving - check if display mode    
    elseif and(GetSecs - checktime > 2.5,GetSecs-triggertime>22) % 2 or 4 s
        checktime = GetSecs;
        bool_updateScreen = true;
        if bool_disponly % simulation mode
            cur_vol = cur_vol+1; rtp.pscraw(cur_vol,:) = (rand(1,2)-0.5)*2;
        elseif istraining == 0 % transfer run
            cur_vol = round((checktime-triggertime-20)/2);
            fprintf('DRIN files not arriving but continue anyways: %0.1fs\n',checktime-triggertime); 
        else
            fprintf('!!WARNING!! DRIN files not arriving: %0.1fs\n',checktime-triggertime); 
            bool_updateScreen = false;
        end
    end
    
% -- Update state of current block
    if any(ds.exp.vols.regu(:)==cur_vol), isrest = false; else, isrest = true; end

% ----------- Display part ------------------------      
    if bool_updateScreen
        fprintf('%.3f vol %d - ',checktime-triggertime,cur_vol);
    % -- rest block: check if baseline needs updating
        if isrest
            if any(ds.exp.vols.rest(end,:)==cur_vol)
                bl = mean(rtp.raw(cur_vol-4:cur_vol,:)); rtp.allbl = [rtp.allbl; bl];
                fprintf('Update baseline\t1: %0.0f, 2: %0.0f', bl(1), bl(2))
            else
                fprintf('Rest block BOLD\t1: %0.0f, 2: %0.0f',rtp.raw(cur_vol,1),rtp.raw(cur_vol,2))
            end
    % -- regulate block: update ptb screen 
        else 
            psc = 100 * (mean(rtp.raw(cur_vol-2:cur_vol,:))-bl)./bl; % if any(isnan(psc)),psc(isnan(psc)) = 0; end
            rtp.pscraw(cur_vol,:) = psc;
            
            if ds.group.db.idx == 1, rtp.dpsc(cur_vol) = rtp.pscraw(cur_vol,1)-rtp.pscraw(cur_vol,2); % [left, right]
            else, rtp.dpsc(cur_vol) = rtp.pscraw(cur_vol,2)-rtp.pscraw(cur_vol,1); end
            fprintf('Regulate BOLD\t1: %0.0f, 2: %0.0f, dPSC: %.2f',rtp.raw(cur_vol,:),rtp.dpsc(cur_vol))
        end
        
    % -- Draw display
        rtp.drawlvl{cur_vol} = draw_NF(pt,isrest,istraining,rtp.dpsc(cur_vol));        
        [VTs, SOT, FTs] = Screen('Flip', pt.win.w);
        ds.r(run).screenflips(end+1).imgArray = Screen('GetImage',pt.win.w); 
        ds.r(run).screenflips(end).timings = [GetSecs, VTs, SOT, FTs];
        if isrest, ds.r(run).screenflips(end).type = sprintf('Rest vol %d',cur_vol);
        else, ds.r(run).screenflips(end).type = sprintf('Regu vol %d',cur_vol);end
    end
end

% -------------------------------------------------------
% Finish scan, clean-up, save data
fprintf('\n%s\tEnd of run %d, runtime: %.1f s\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),run,GetSecs-triggertime);
if eye
    try
        NF_eyelink('stopr',edfFile); ds.r(run).el = el; 
    catch
        disp('Did not close edfFile or stopped eye tracking')
    end
end
Screen('closeall');

ds.r(run).isdone = 1; ds.r(run).pt = pt; ds.r(run).rtp = rtp;
ds.r(run).screenflips = struct2table(ds.r(run).screenflips);
fprintf('\n\n----------------------\nCheck if DDD.exe has finished transferring all files (vols 3, 4, 5)\n\n')
ds.r(run).comments = input('Comments about this run: ','s');
if ~bool_disponly, anonymize_files(mr8,ds,run);end

% -- Prep next run or finish session
fin = input('Finished the session? y/[n]: ','s');
if strcmp(fin,'y')==0
    fprintf('Check motion in TBV\nReady for next run!-> run current section\n')
else
    fprintf('Final run, session finished!\nSaving ds.mat file\n')
    if exist([ds.dirs.ps filesep 'ds.mat'],'file')==0
        save([ds.dirs.ps filesep 'ds.mat'],'ds','-v7.3');
        disp('ds.mat saved')
    else
        fprintf('%s already exists!\nSaved as ds%d.mat\n',[ds.dirs.ps filesep 'ds.mat'],length(dir([ds.dirs.ps filesep 'd*.mat']))+1)
        save([ds.dirs.ps filesep 'ds' num2str(length(dir([ds.dirs.ps filesep 'd*.mat']))+1) '.mat'],'ds','-v7.3');
    end
    if eye, Eyelink('Shutdown'); end
    diary off
%     Priority(0); % return to normal priority
end

%% After transferring EDF files through popup_calibration
disp('Move to the folder with the EDF files')
edfFiles = dir('*.edf');

try
    for i = 1:length(edfFiles),system(sprintf('edf2asc %s',edfFiles(i).name)); end
    ascFiles = dir('*.asc'); for i = 1:length(ascFiles),movefile(ascFiles(i).name,ds.dirs.ps);end
catch
    disp('Error in conversion edf2asc, edfFiles not moved.')
end
for i = 1:length(edfFiles),movefile(edfFiles(i).name,ds.dirs.ps);end

    
    