% PROJECT:      WP3 - rt-fMRI NF for self-regulation of interhemispheric IPS activity
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Processing of all NIFTI data across all participants
% 1. Anatomical: segment structural file of first NF session
%                reslice to 1mm (visualisation), 2.5mm (task-related), 4mm
%                (functional connectivity): res1_nii, res2_nii, res4_nii
% 2. Functional: re-align, coregister to 3D, smooth
%                check outlier with ART toolbox (manual)
% 4. FC:         pre-process resting-state volumes
% OUTPUT:       Processed NIFTI data in /DataDerived + ART outliers
% -------------------------------------------------------------------------
% 2020.06.04 Stay in native space (v1_1)
% 2020.06.05 Fix co-registration (v1_2)
% 2020.06.25 Change to remove first 10 volumes (v2)

function dirs = WP3_2_preproc_nii(dirs,qc)
if nargin<2, qc=false; end

% Initialisation
fprintf('\n\n***************** Start WP3_2_preproc_nii *****************\n\n')

%% 1. Anatomical 

% 1. Get data
fprintf('Structural files: %d/%d\n',numel(dirs.data.anat),dirs.n.p*dirs.n.s)

% 2. Check segmentation
if 0
for p = 1:dirs.n.p
    spmdir = spm('Dir');
    structfol = dirs.data.anat{s,p};
    if exist([structfol '/iy_3D.nii'],'file')==0
        fprintf('\nPerforming segmentation on %s\n',structfol); 
        clear matlabbatch
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {[structfol '3D.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spmdir '/tpm/TPM.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spmdir '/tpm/TPM.nii,2']};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spmdir '/tpm/TPM.nii,3']};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spmdir '/tpm/TPM.nii,4']};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spmdir '/tpm/TPM.nii,5']};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spmdir '/tpm/TPM.nii,6']};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 1; 
        matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
        spm_jobman('run',matlabbatch); 
    end
    fprintf('Finished segmentation of %s!\n',structfol);
    % Normalise to several resolutions
    if p == 1, vres = [1 2.5 4];else, vres = [2.5 4];end
    for i = 1:length(vres)
        if exist([structfol '/w' num2str(floor(vres(i))) '_3D.nii'],'file')==0
            fprintf('Performing normalisation of %s...\n',structfol)
            clear matlabbatch
            matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[structfol 'y_3D.nii']};
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[structfol '3D.nii']};
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                      78 76 85];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = ones(1,3)*vres(i);
            matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
            spm_jobman('run',matlabbatch);

            movefile([structfol '/w3D.nii'],[structfol 'w' num2str(floor(vres(i))) '_3D.nii'])
        end
        disp(['Finished reslicing ' structfol 'w' num2str(floor(vres(i))) '_3D.nii'])
    end
    fprintf('%s segmentation complete!\n\n',dirs.data.anat{p}(length(dirs.data.main):end))
end
end
disp(' ')


%% 2. Functional 

% 1. Get data
fprintf('Functional files: %d/%d\n\n',numel(dirs.data.func),dirs.n.p*dirs.n.r*dirs.n.s)


% 2. Check pre-processing
dirs.processed.func = cell(size(dirs.data.func));
dirs.processed.funcfn = ['sr' dirs.data.funcfn];

spm('defaults', 'FMRI');
v_new = 2.5;

for p = 1%:dirs.n.p
    for s = 1:dirs.n.s
        for r = 1:dirs.n.r
            fn_result = [dirs.data.func{r,s,p} dirs.processed.funcfn];
            
            if exist(fn_result,'file')==0
                % pre-processing task-related functional files
                ff = [dirs.data.func{r,s,p} dirs.data.funcfn]; 
                [funcfol,fname,ext] = fileparts(ff);
                structfol = dirs.data.anat{s,p};
                
                % 4D 220->210 volumes; removing first 10 volumes, keep the header but 
                % change the orientation in V(i).mat
                skip_images=10;
                system(sprintf('chflags -R nouchg %s',ff));
                V=spm_vol(ff); nvol=length(V);
                if nvol == 220
                fprintf('Skip 10 volumes %s...\n',ff)
                data=zeros(V(1).dim(1),V(1).dim(2),V(1).dim(3),nvol-skip_images);
                for i=1:nvol, data(:,:,:,i)=spm_read_vols(V(i)); end
                delete(ff);
                Vt=V(1:nvol-skip_images);
                for i=1:nvol-skip_images
                    Vt(i).mat=V(i+skip_images).mat;
                    spm_write_vol(Vt(i),data(:,:,:,i+skip_images));
                end
                end
                
                fprintf('Start pre-processing %s...\n',ff)

                fprintf('realigning %s...\n',ff) 
                if or(and(r == 1,exist(strrep(ff,'4D','mean4D'),'file')==0),...
                        exist([funcfol '/rp_' fname '.txt'],'file')==0)
                    clear matlabbatch
                    if r == 1
                    matlabbatch{1}.spm.spatial.realign.estwrite.data = {{ff}}';
                    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 2;
                    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
                    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
                    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
                    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1]; % mean only
                    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
                    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
                    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
                    elseif r>1
                    matlabbatch{1}.spm.spatial.realign.estimate.data = {{ff}}';
                    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
                    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 2;
                    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
                    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 1;
                    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
                    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
                    end
                    spm_jobman('run',matlabbatch);
                end

                % Coregister anatomical to mean functional
                if r == 1
                    clear matlabbatch
                    matlabbatch{1}.spm.spatial.coreg.write.ref = {strrep(ff,'4D','mean4D')};
                    matlabbatch{1}.spm.spatial.coreg.write.source = {[structfol '/3D.nii,1']};
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
                    spm_jobman('run',matlabbatch);
                end
                fprintf('coregistering %s...\n',ff)
                
                if exist([structfol '/r4D.nii'],'file')==0
                    clear matlabbatch
                    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[structfol '/3D.nii,1']};
                    matlabbatch{1}.spm.spatial.coreg.estimate.source = {[ff ',1']};
                    matlabbatch{1}.spm.spatial.coreg.estimate.other = {ff};
                    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = ...
                        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                    spm_jobman('run',matlabbatch);

                    % COREGISTER write
                    clear matlabbatch
                    matlabbatch{1}.spm.spatial.coreg.write.ref = {[structfol '/r3D.nii,1']};
                    matlabbatch{1}.spm.spatial.coreg.write.source = {ff};
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
                    spm_jobman('run',matlabbatch);
                end

                fprintf('smoothing %s...\n',ff)
                if exist([structfol '/' dirs.processed.funcfn],'file')==0
                    clear matlabbatch
                    matlabbatch{1}.spm.spatial.smooth.data = {[funcfol '/r' fname ext]};
                    matlabbatch{1}.spm.spatial.smooth.fwhm = 2*v_new*ones(1,3);
                    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                    matlabbatch{1}.spm.spatial.smooth.im = 0;
                    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

                    spm_jobman('run',matlabbatch);
                end
            end
            fprintf('%s preprocessing complete!\n',fn_result(length(dirs.data.main):end))
            dirs.processed.func{r,s,p} = fn_result;
        end
    end
    disp(' ')
end

%% qc: check max motion
fprintf('\nCheck maximum head motion\n')
for p = 1%:dirs.n.p
    fprintf('\nParticipant %d\n',p)
    for s = 1:dirs.n.s
        fprintf('Session %d\n',s)
        for r = 1:dirs.n.r
            R = load([fileparts(dirs.processed.func{r,s,p}) '/rp_4D.txt']);
            R(:,4:6) = rad2deg(R(:,4:6)); R = R(:,1:6); 
            fprintf('r%d max: %.2f, intervol: %.2f',...
                r,max(abs(R(:))),max(max(abs(diff(R)))))
            if max(abs(R(:)))>2.5
                fprintf('\t!!!\n')
            else
                fprintf('\n')
            end
        end
    end
end


%% 3. Check outliers with ART

if 1%qc
    dirs.data.art = reshape(cellfun(@(x) strrep(x,'run-1','art'),...
        dirs.data.func(1,:,:),'UniformOutput',false),dirs.n.s,dirs.n.p);
dirs.data.artreg = cell(5,3,6);
    for p = 1%:dirs.n.p
        for s = 1:dirs.n.s
            sess_file = [dirs.data.art{s,p} '/artsettings.cfg']; 
            if exist(sess_file,'file')==0
                mkdir(dirs.data.art{s,p});cd(dirs.data.art{s,p})
                % Write settings file (edited from art_batch.m)
                fid=fopen(sess_file,'wt');
                fprintf(fid,'sessions: %d\n',dirs.n.r);
                fprintf(fid,'global_mean: %d\n',1);
                fprintf(fid,'global_threshold: %f\n',5.0);
                fprintf(fid,'motion_threshold: %f\n',2.0);
                fprintf(fid,'motion_file_type: %d\n',0);
                fprintf(fid,'motion_fname_from_image_fname: 1\nend\n');
                for n2=1:dirs.n.r, fprintf(fid,'session %d image %s4D.nii\n',n2,dirs.data.func{n2,s,p}); end
                fprintf(fid,'end\n');
                fclose(fid);
            end
            fn_artreg = dir(sprintf('%srun-*/art_regression_outliers_and_movement_4D.mat',...
                dirs.data.art{s,p}(1:end-4)));
            if isempty(fn_artreg)
                % perform ART outlier detection
                fprintf('\nPerform %s ART outlier detection\n',...
                    dirs.data.art{s,p}(length(dirs.data.main):end));tic
                art('sess_file',sess_file);t=toc;
                fprintf('Elapsed time: %dm%ds\n',floor(t/60),round(mod(t,60)));
                pause; close all;
                fn_artreg = dir(sprintf('%srun-*/art_regression_outliers_and_movement_4D.mat',...
                    dirs.data.art{s,p}(1:end-4)));
            end
            fprintf('%s ART outlier detection complete!\n',...
                dirs.data.art{s,p}(length(dirs.data.main):end))
            dirs.data.artreg(:,s,p) = fullfile({fn_artreg.folder}',{fn_artreg.name}');
        end
    end
end

%% End 
fprintf('\n\n***************** End WP3_2_preproc_nii *****************\n\n')



