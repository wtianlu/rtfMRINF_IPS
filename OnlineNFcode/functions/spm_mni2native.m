% 2019.09.23 Deformations: set image defining fov to struct3D instead of y_struct
% 2019.07.30 matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm from 40 to 60
%            add moving files here
% 2019.06.26 With new MATLAB environment
% 2019.06.21 In new MR8 & Python environment
% 2018.03.06 removed conversion to Turbo Brainvoyager format to separate .m
% 2017.11.08 set threshold to 0.5; change name to include 'rw'
% 2017.08.31 first version

%% Select and move file
function ds = spm_mni2native(ds)
fprintf('***spm_mni2native.m************************\n')
spm_dir = spm('Dir');
spm('defaults','fMRI');

if ds.idstr(end)=='1'
    nv = length(dir([ds.dirs.struct filesep 'voi*.nii']));

    if exist([ds.dirs.struct filesep 'wvoi1.nii'],'file')==0
        disp('Start segmentation of struct3D.nii...')
        clear matlabbatch;
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {[ds.dirs.struct filesep 'struct3D.nii']};
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
        spm_jobman('run',matlabbatch);
        disp('Finished struct3D.nii segmentation!')

        disp('Start warping template VOIs to native space...')
        fp_vois = cellfun(@(x) [ds.dirs.struct filesep 'voi' num2str(x) '.nii'],num2cell(1:nv),'UniformOutput',false)';

        clear matlabbatch
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[ds.dirs.struct filesep 'y_struct3D.nii']}; % push y or pull iy -> interpolation or not
        matlabbatch{1}.spm.util.defs.out{1}.push.fnames = fp_vois;
        matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
        matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {ds.dirs.struct};
        matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {[ds.dirs.struct filesep 'struct3D.nii']}; % image defining file same as struct
        matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
        matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.push.prefix = 'w';
        spm_jobman('run',matlabbatch);
        disp('Finished warping template VOIs to native space!')
    else
        fprintf('%s already exists!\n',[ds.dirs.struct filesep 'wvoi1.nii']);
    end
else
    disp('Select structural from session 1')
    [struct_s1,struct_s1dir] = uigetfile('*.nii');disp(fullfile(struct_s1dir,struct_s1));
    copyfile(fullfile(struct_s1dir,struct_s1), [ds.dirs.struct filesep 's1_struct3D.nii'])
    for i = 1:2
        fprintf('Select ROI %d from previous session\n',i)
        [rois_s1{i},rois_s1dir] = uigetfile('*.nii');
        copyfile(fullfile(rois_s1dir,rois_s1{i}),[ds.dirs.struct filesep 's1_' rois_s1{i}])
    end
    
    clear matlabbatch
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ds.dirs.struct filesep 'struct3D.nii']}; 
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[ds.dirs.struct filesep 's1_struct3D.nii']}; 
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellfun(@(x) [ds.dirs.struct filesep 's1_' x],rois_s1,'UniformOutput',false)'; % let VOIs follow
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);

    disp('Finished coregistration to current session!')
    for i = 1:2
        movefile([ds.dirs.struct filesep 'rs1_brwvoi' num2str(i) '.nii'],[ds.dirs.struct filesep 'wvoi' num2str(i) '.nii'])
    end
    movefile([ds.dirs.struct filesep 'rs1_struct3D.nii'],[ds.dirs.struct filesep 'rstruct3D.nii'])
    nv = length(dir([ds.dirs.struct filesep 'wvoi*.nii']));
end
    
% Coregister to functional space    
if exist([ds.dirs.rois filesep 'brwvoi1.nii'],'file')==0
    disp('Start coregistration VOIs to functional space...')
    fp_roinative = cellfun(@(x) [ds.dirs.struct filesep 'wvoi' num2str(x) '.nii,1'],num2cell(1:nv),'UniformOutput',false)';

    clear matlabbatch
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ds.dirs.rois filesep 'func4D_check.nii,1']}; % functional file
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[ds.dirs.struct filesep 'struct3D.nii,1']}; % transform according to higher res struct
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = fp_roinative; % let VOIs follow
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;%(str2double(ds.idstr(end))>1)*4; % 0 if s1, 4 if s2, s3
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);

    disp('Finished coregistration to functional space!')

    % binarize coregistered vois
    clear matlabbatch
    fp_roi = cellfun(@(x) [ds.dirs.struct filesep 'rwvoi' num2str(x) '.nii,1'],num2cell(1:nv),'UniformOutput',false)';
    for v = 1:length(fp_roi)
        matlabbatch{v}.spm.util.imcalc.input = fp_roi(v);
        matlabbatch{v}.spm.util.imcalc.output = [ds.dirs.struct filesep 'brwvoi' num2str(v) '.nii'];
        matlabbatch{v}.spm.util.imcalc.outdir = {''};
        matlabbatch{v}.spm.util.imcalc.expression = 'i1>0.5';
        matlabbatch{v}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{v}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{v}.spm.util.imcalc.options.mask = 0;
        matlabbatch{v}.spm.util.imcalc.options.interp = 1;
        matlabbatch{v}.spm.util.imcalc.options.dtype = 4;
    end
    spm_jobman('run',matlabbatch);

    movefile([ds.dirs.struct filesep 'rstruct3D.nii'],[ds.dirs.rois filesep 'rstruct3D.nii'])
    for v = 1:nv
        movefile([ds.dirs.struct filesep 'rwvoi' num2str(v) '.nii'],[ds.dirs.rois filesep 'rwvoi' num2str(v) '.nii'])
        movefile([ds.dirs.struct filesep 'brwvoi' num2str(v) '.nii'],[ds.dirs.rois filesep 'brwvoi' num2str(v) '.nii'])
    end
else
    fprintf('%s already exists!\n',[ds.dirs.rois filesep 'brwvoi1.nii']);
end

fprintf('*******************************************\n')











