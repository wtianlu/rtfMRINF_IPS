% PROJECT:      WP3 - rt-fMRI NF for self-regulation of interhemispheric IPS activity
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Prepare all NIFTI data across all participants for pre-processing
%   1. Raw data:    check complete nii dataset: 3D, 4D rest and task
%   2. QC:          view volumes, check for missing data
%   3. Prepare:     copy raw files to perform pre-processing
%                   get conditions and online data
% OUTPUT:       NIFTI data organisedin /DataDerived
% -------------------------------------------------------------------------
% 2020.06.25 Correct conditions.mat to include onset+durations of rest
% blocks in order to disentangle effect from initial fixation period
% 2020.06.05 Change to include anatomical scan from every session (v1_2)
function dirs = WP3_1_prepfiles_nii(dirs,qc)

% Initialisation
if nargin<2, qc=false; end
fprintf('\n\n***************** Start WP3_1_prepfiles_nii *****************\n\n')

% prepare filenames
dirs.data.anatfn = '3D.nii';
dirs.data.funcfn = '4D.nii';
dirs.data.restfn = '4D.nii';
dirs.data.condfn = 'conditions.mat';

%% Check anatomical data
fprintf('Check anatomical data\n')
anat = dir([dirs.raw.main 'sub-*/ses-*/anat/subj*_s*_T1w.nii']);
if length(anat)~=dirs.n.s*dirs.n.p,disp('Missing T1w data!');disp(anat);pause;end
anat = reshape(fullfile({anat.folder}',{anat.name}'),dirs.n.s,dirs.n.p);
dirs.raw.anat = anat;

if qc
fprintf('Quality check structural data...\n\n')
for p = 1:dirs.n.p
    spm_check_registration(anat{:,p})
    pause; close all;
end
end

%% Check functional data
fprintf('Check functional data\n')
func = dir([dirs.raw.main 'sub-*/ses-*/func/subj*_s*_task-*.nii']);
if length(func)~=dirs.n.r*dirs.n.s*dirs.n.p,disp('Missing task data!');disp(func);pause;end
func = reshape(fullfile({func.folder}',{func.name}'),dirs.n.r,dirs.n.s,dirs.n.p);
dirs.raw.func = func;

if qc
fprintf('Quality check functional data...\n\n')
for p = 1:dirs.n.p
    for s = 1:dirs.n.s
        tmp = cellfun(@(x) [x ',1'],func(:,s,p),'UniformOutput',false);
        spm_check_registration(tmp{:})
        pause; close all;
        
        for r = 1:dirs.n.r
            if ~isdeployed, addpath(fullfile(spm('Dir'),'spm_orthviews')); end
            spm_ov_browser('ui',spm_select('expand',func{r,s,p}));
            pause; close all
        end
        
    end
end
end

%% Check resting-state data
fprintf('Check resting-state data\n')
rest = dir([dirs.raw.main 'sub-*/ses-*/func/subj*_s*_rest.nii']);
if length(rest)~=2*dirs.n.p,disp('Missing rest data!');disp(rest);pause;end
rest = reshape(fullfile({rest.folder}',{rest.name}'),dirs.n.rest,dirs.n.p);

dirs.raw.rest = rest;

if qc
fprintf('Quality check resting-state data...\n\n')
for p = 1:dirs.n.p
    tmp = cellfun(@(x) [x ',1'],rest(:,p),'UniformOutput',false);
    spm_check_registration(tmp{:})
    pause; close all;
        
    for s = 1:size(rest,1)
        if ~isdeployed, addpath(fullfile(spm('Dir'),'spm_orthviews')); end
        spm_ov_browser('ui',spm_select('expand',rest{s,p}));
        pause; close all

    end
end
end

%% Copy files to processing directory

fprintf('Prepare NIFTI data for pre-processing...\n\n')

% define directories
dirs.data.main = [dirs.main 'DataDerived/'];
dirs.data.anat = cell(dirs.n.s,dirs.n.p);
dirs.data.func = cell(dirs.n.r,dirs.n.s,dirs.n.p);
dirs.data.rest = cell(dirs.n.rest,dirs.n.p);


for p = 1:dirs.n.p
    
    for s = 1:dirs.n.s
        % Structural nii
        dirs.data.anat{s,p} = [dirs.data.main 'sub-0' num2str(p) '/anat/ses-' num2str(s) '/'];
        if exist(dirs.data.anat{s,p},'dir')==0,mkdir(dirs.data.anat{s,p});end
        if exist([dirs.data.anat{s,p} dirs.data.anatfn],'file')==0
            copyfile(dirs.raw.anat{1,p},[dirs.data.anat{s,p} dirs.data.anatfn]);
        end
        fprintf(['\nSaved ' dirs.data.anat{s,p} dirs.data.anatfn '\n'])
    
        % Functional nii
        for r = 1:dirs.n.r
            dirs.data.func{r,s,p} = sprintf('%ssub-0%d/func/ses-%d/run-%d/',dirs.data.main,p,s,r);
            if exist(dirs.data.func{r,s,p},'dir')==0,mkdir(dirs.data.func{r,s,p});end
            if exist([dirs.data.func{r,s,p} dirs.data.funcfn],'file')==0
                copyfile(dirs.raw.func{r,s,p},[dirs.data.func{r,s,p} dirs.data.funcfn]);
            end
            disp(['Saved ' dirs.data.func{r,s,p} dirs.data.funcfn])
        end
        
        % Rest nii
        if s <= dirs.n.rest
            dirs.data.rest{s,p} = sprintf('%ssub-0%d/rest/ses-%d/',dirs.data.main,p,s);
            if exist(dirs.data.rest{s,p},'dir')==0,mkdir(dirs.data.rest{s,p});end
            if exist([dirs.data.rest{s,p} dirs.data.restfn],'file')==0
                copyfile(dirs.raw.rest{s,p},[dirs.data.rest{s,p} dirs.data.restfn]);
            end
            disp(['Saved ' dirs.data.rest{s,p} dirs.data.restfn])
        end
    end
end

%% Read diary data and get conditions and dPSC

fn_online = [dirs.results.main 'NF_online_dPSC.mat'];

if exist(fn_online,'file')>0
    load(fn_online);disp('online.mat loaded!')
else
    fprintf('\nCreate online.mat\n')
    dirs.online.diary = dir([dirs.raw.main 'sub-*/ses-*/online/NSL*_diary']);
    dirs.online.diary = fullfile({dirs.online.diary.folder}',{dirs.online.diary.name}');
    for p = 1:dirs.n.p
        for s = 1:dirs.n.s
            % --- read diary file
            fndiary = dirs.online.diary{contains(dirs.online.diary,sprintf('sub-0%d/ses-%d',p,s+1))};
            disp(fndiary)
            fid = fopen(fndiary,'r'); 
            tline = fgetl(fid);tlines = {};
            while ischar(tline)
                tlines{end+1} = tline;
                tline = fgetl(fid);
            end
            fclose(fid); tlines = tlines';
            clear fid tline;

            % --- extract data
            idx = [find(contains(tlines,'Trigger arrived for run')); length(tlines)];
            run_diary = cellfun(@(x1,x2) tlines(x1:x2), num2cell(idx(1:end-1)),...
                num2cell(idx(2:end)),'UniformOutput',false);
            for r = 1:length(run_diary)%nr
                vol_data = cellfun(@strsplit,run_diary{r}(contains(run_diary{r},'- level:')),...
                    'UniformOutput',false);
                vol_nr = cellfun(@(x) str2double(x(3)),vol_data);
                vol_regu = find(cellfun(@(x) ~isempty(find(contains(x,'dPSC'),1)),vol_data));
                screenflips = nan(210,1);
                screenflips(vol_nr) = str2double(cellfun(@(x) x(1),vol_data));
                tlevel = nan(210,1);
                tmp = cellfun(@(x) x(contains(x,'/13')),vol_data);
                tmp = cellfun(@(x) strsplit(x,'/'),tmp,'UniformOutput',false);
                tlevel(vol_nr) = str2double(cellfun(@(x) x(1),tmp));
                bold = nan(210,2);
                bold(vol_nr,1) = cellfun(@(y) str2double(y(1:end-1)),...
                    cellfun(@(x) x(find(strcmp(x,'1:'))+1),vol_data));
                tmp = cellfun(@(x) x(find(strcmp(x,'2:'))+1),vol_data);
                tmp(cellfun(@(x) contains(x,','),tmp)) = cellfun(@(x) x(1:end-1),...
                    tmp(cellfun(@(x) contains(x,','),tmp)),'UniformOutput',false);
                bold(vol_nr,2) = str2double(tmp);
                dpsc = nan(210,1);
                dpsc(vol_regu) = cellfun(@(x) str2double(x(12)),vol_data(vol_regu));

                online.diarysf{r,s,p} = screenflips;
                online.diarytlevel{r,s,p} = tlevel;
                online.diarybold{r,s,p} = bold; % not swapped: 1: left, 2: right
                online.diarydpsc{r,s,p} = dpsc;
                
                % --- Save conditions.mat
                fncond = [dirs.data.func{r,s,p} dirs.data.condfn];
                if 1%exist(fncond,'file')==0
                onsets = screenflips(11:25:end)-20;
                if any(isnan(onsets)), onsets(isnan(onsets)) = ...
                            onsets(find(isnan(onsets))-1)+50;end;
                durations = screenflips(26:25:end)-20-onsets;
                if any(isnan(durations)), durations(isnan(durations)) = ...
                            durations(find(isnan(durations))-1);end
                names = {'regulate'};
                save(fncond, 'names','onsets','durations')
                end
                fprintf('Saved %s\n',[dirs.data.func{r,s,p} dirs.data.condfn])
            end
            online.diaryraw{:,s,p} = run_diary;

        end
    end
    save(fn_online,'online');
    fprintf('\nSaved %s!\n',fn_online)
end


fprintf('\n\n***************** End WP3_1_prepfiles_nii *****************\n\n')








