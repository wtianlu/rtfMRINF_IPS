
function ds = initialise_nf(mr8,whichstep,ds,idP,idS,run)
fprintf('***initialise_nf.m %d***********************\n',whichstep)
bool_testing = ~mr8; %if not at mr8, then I'm testing and copying files, else move out of export folder

switch whichstep
    case 0
        fprintf('%s\tInitializing participant %s session %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),idP, idS)

        % get paths and directories
        dirs.code = pwd;
        [dirs.main,~,~] = fileparts(dirs.code);
        dirs.tpl = [dirs.code filesep 'Templates'];
        dirs.p = [dirs.main filesep 'Data' filesep 'NSL' idP];  
        dirs.ps = [dirs.p filesep 's' idS]; 
        if exist(dirs.ps,'dir')~=0 
            fprintf('\n%s already exists!\nPress key to continue with current idP/idS or cancel code execution\n',dirs.ps); 
            pause;
        end
        dirs.struct = [dirs.ps filesep 'struct']; % structural data in session 1 only
        dirs.tbvs = [dirs.ps filesep 'tbvsettings']; 
        dirs.rois = [dirs.ps filesep 'roi']; 
        
        if exist(dirs.ps,'dir')==0
            mkdir(dirs.ps);mkdir(dirs.struct); mkdir(dirs.tbvs);
            mkdir(dirs.rois)
            disp('Participant session directories created!')
        else
            disp('Participant session directory already exists'); pause
        end 
        
        % participant and experiment settings
        ds.idstr = sprintf('NSL%sS%s',idP,idS); 
        ds.exp.drinfirstfname = '*.img';
        ds.exp.nvol = 220; ds.exp.nrun = 5;
        allvols = reshape(1:(ds.exp.nvol)+5,25,numel(1:(ds.exp.nvol)+5)/25);
        ds.exp.vols = struct('rest',allvols(1:10,:),'regu', allvols(11:end,1:end-1));

        % Prepare ROI templates in session 1
        if idS == '1'
            if exist([dirs.struct filesep 'voi1.nii'],'file')==0
                tpl_voi = dir([dirs.tpl filesep 'VOI' filesep 'tpl_*.nii']); 
                for vv = 1:length(tpl_voi)
                    copyfile([dirs.tpl filesep 'VOI'  filesep tpl_voi(vv).name],[dirs.struct filesep 'voi' num2str(vv) '.nii']); 
                    fprintf('%d: voi%s\n', vv,tpl_voi(vv).name(4:end))
                end
                fprintf('Copied %d voi template files\n',vv)
            else
                disp('VOI template files already present')
            end
        end
        ds.dirs = dirs;
% ------------------------------------------------------------------------
    case 1
        dirs = ds.dirs;
        % select neuroimaging files
        if exist([dirs.struct '/struct3D.nii'],'file')==0
            disp('Select structural file...');
            [fp_struct, dirs.structtmp] = uigetfile('*.nii', 'Select the structural file');
            if fp_struct~=0
                if bool_testing, copyfile([dirs.structtmp filesep fp_struct],[dirs.struct filesep 'struct3D.nii']);
                else,            movefile([dirs.structtmp filesep fp_struct],[dirs.struct filesep 'struct3D.nii']);
                end
                fprintf('%s\t->\t%s\n',[dirs.structtmp filesep fp_struct],[dirs.struct filesep 'struct3D.nii']);
            else
                fprintf('No structural file selected!\n'); return;
            end                
        else
            fprintf('%s already present!\n',[ds.dirs.struct '/struct3D.nii']);
        end
        
        if exist([dirs.rois '/func4D_check.nii'],'file')==0
            disp('Select functional_check file...'); 
            [fp_func, dirs.structtmp] = uigetfile('*.nii', 'Select the functional file');
            if fp_func~=0
                if bool_testing, copyfile([dirs.structtmp filesep fp_func],[dirs.rois filesep 'func4D_check.nii']);
                else,            movefile([dirs.structtmp filesep fp_func],[dirs.rois filesep 'func4D_check.nii']);
                end
                fprintf('%s\t->\t%s\n',[dirs.structtmp filesep fp_func],[dirs.rois filesep 'func4D_check.nii'])
            else
                fprintf('No Functional file selected!\n')
            end
        else
            fprintf('%s already present!\n',[ds.dirs.rois '/func4D_check.nii'])
        end
        if isempty(dir([dirs.rois '/*.img']))
            disp('Select sample DRIN.img file from functional_check...'); 
            [fn_pdrin,fpdrin] = uigetfile('*.img', 'Select the DRIN file'); 
            if bool_testing
                copyfile([fpdrin filesep fn_pdrin],[dirs.rois filesep fn_pdrin]);
            else
                movefile([fpdrin filesep fn_pdrin],[dirs.rois filesep fn_pdrin]);
            end
            fprintf('Move %s from %s\t->\t%s\n',fn_pdrin,fpdrin,dirs.rois)
            fn_pdrin = strsplit(fn_pdrin,{'-','.'});
            ds.exp.drinfirstfname = fn_pdrin{1};
            disp('Finished moving the brain files!')
        else
            pdrin = dir([dirs.rois '/*.img']);
            fprintf('%s already present!\n',pdrin.name)
            fn_pdrin = strsplit(pdrin.name,{'-','.'});
            ds.exp.drinfirstfname = fn_pdrin{1};
        end
        
% ------------------------------------------------------------------------
    case 2
        dirs = ds.dirs;
        fprintf('Setting up directories for run %d...\n',run);
        ds.r(run).dir = [dirs.ps filesep 'r' num2str(run)];
        ds.r(run).timings = struct('mrtrigger',[],'rtp', [],'drin',[]);
        ds.r(run).rtp = struct('raw',zeros((ds.exp.nvol)-10,2), 'pscraw', zeros((ds.exp.nvol)-10,2),...
            'dpsc',zeros((ds.exp.nvol)-10,1),'drawlvl', [], 'allbl', []);
        ds.r(run).screenflips = struct('imgArray',{},'timings',{},'type',{});
        ds.r(run).isdone = 0;
        if exist(ds.r(run).dir,'dir')==0
            mkdir(ds.r(run).dir);
            mkdir([ds.r(run).dir filesep 'tbv_feedback']);
            mkdir([ds.r(run).dir filesep 'tbv_target']);
            mkdir([ds.r(run).dir filesep 'tbv_watch']);
        end

        fprintf('Writing TBV settings files for run %d...\n',run);
        fp_tbvprt = [dirs.tpl filesep 'TBV' filesep 'WP3p_continuous_220.prt'];
        disp('Select new DRIN folder')
        fp_tbvwatchfolder = uigetdir('Z:/datadump/');
        ds.dirs.exportdrin = fp_tbvwatchfolder;
        fp_tbvsettings = [dirs.tpl filesep 'TBV' filesep 'WP3p_template.tbv'];
        tbv_template_str = fileread(fp_tbvsettings);

        tbv_data = {sprintf('%s_r%d',ds.idstr,run);     % Title
            fp_tbvwatchfolder;                          % watchfolder
            [ds.r(run).dir filesep 'tbv_target'];       % targetfolder
%             ds.exp.drinfirstfname;                      % firstfilename 
            num2str(ds.exp.nvol);                       % number of total volumes
            [ds.idstr '_r' num2str(run)];               % slice prefix
            [ds.idstr '_r' num2str(run)];               % project name
            fp_tbvprt;                                  % stimulation protocol 
            [ds.r(run).dir filesep 'tbv_feedback']      % feedback folder
            };
        fid = fopen([dirs.tbvs filesep ds.idstr '-' num2str(run) '.tbv'],'w');
        fprintf(fid,tbv_template_str,tbv_data{:});
        fclose(fid);  

        disp('Finished writing TBV settings file!')
end
fprintf('*******************************************\n')






