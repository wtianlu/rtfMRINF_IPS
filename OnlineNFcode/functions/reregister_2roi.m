% Re-coregister ROIs if participant moved between runs and transform to
% *.roi format in TBV

%% Get files
disp('Select new fMRI_check4.nii')
[fmricheck_new,fmricheck_dir] = uigetfile('*.nii');disp(fullfile(fmricheck_dir,fmricheck_new));

% outdir = [ds.dirs.rois num2str(length(dir([ds.dirs.rois '*']))+1)];
outdir = pwd;
if exist(outdir,'dir')==0,mkdir(outdir);end
fprintf('Move new fMRI_check4.nii to %s\n',outdir)
movefile(fullfile(fmricheck_dir,fmricheck_new),[outdir filesep 'func4D_check.nii'])
fmricheck_new = 'func4D_check.nii';

disp('Select old fMRI_check4.nii')
[fmricheck_old,fmricheck_dir] = uigetfile('*.nii');disp(fullfile(fmricheck_dir,fmricheck_old));
fprintf('Copy old fMRI_check4.nii to %s\n',outdir)
copyfile(fullfile(fmricheck_dir,fmricheck_old),[outdir filesep 'func4D_check_old.nii'])
fmricheck_old = 'func4D_check_old.nii';

rois_s1 = cell(1,2);
for i = 1:2
    fprintf('Select old ROI %d from current session\n',i)
    [rois_s1{i},rois_s1dir] = uigetfile('*.nii');
    copyfile(fullfile(rois_s1dir,rois_s1{i}),[outdir filesep 'old_' rois_s1{i}])
end
disp('Finished prepping files')

%% Coregister previous fMRI_check4 to new one, let ROIs move along

clear matlabbatch
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr([outdir filesep fmricheck_new ',1']); 
matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr([outdir filesep fmricheck_old ',1']); 
matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(cellfun(@(x) [outdir filesep 'old_' x],rois_s1,'UniformOutput',false)'); % let VOIs follow
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
disp('Finished coregistration to new fMRI4_check!')
    
% binarize coregistered vois
clear matlabbatch
matlabbatch = cell(1,2);
for v = 1:2
    matlabbatch{v}.spm.util.imcalc.input = {[outdir filesep 'rold_' rois_s1{v}]};
    matlabbatch{v}.spm.util.imcalc.output = ['brwvoi' num2str(v) '.nii'];
    matlabbatch{v}.spm.util.imcalc.outdir = {outdir};
    matlabbatch{v}.spm.util.imcalc.expression = 'i1>0.5';
    matlabbatch{v}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{v}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{v}.spm.util.imcalc.options.mask = 0;
    matlabbatch{v}.spm.util.imcalc.options.interp = 1;
    matlabbatch{v}.spm.util.imcalc.options.dtype = 4;
end
spm_jobman('run',matlabbatch);

disp('Finished binarising new brwvois!')

%% Convert to .roi format

% Copy template 1 
fn_result = [outdir filesep ds.idstr '_tbv_voi' outdir(end) '.roi'];
copyfile([ds.dirs.tpl filesep 'TBV' filesep 'ROIs_template1.roi'],fn_result);

fp_roi = cellfun(@(x) [outdir filesep 'brwvoi' num2str(x) '.nii'],num2cell(1:2),'UniformOutput',false)';
for ii = 1:length(fp_roi)
    fprintf('Loading %s\n',fp_roi{ii})
    % load nifti
    V = spm_vol(fp_roi{ii}); [dv,~] = spm_read_vols(V);
    % get coordinates
    ind=find(dv>0.8); [i1, i2, i3] = ind2sub(size(dv), ind); 
    xyz{ii} = [i1, i2, i3]'; nVoxels = length(i1);
    fprintf(' nVoxels: %d x: %d-%d y: %d-%d slices: %d-%d\n',nVoxels,min(i1),max(i1),min(i2),max(i2),min(i3),max(i3))
    
    % Get header info
    Zslice = round(mean(i3));
    Xleft = min(i1(i3==Zslice)); Xright = max(i1(i3==Zslice)); 
    Ytop = min(i2(i3==Zslice)); Ybottom = max(i2(i3==Zslice)); 
    ind3D = [i1, i2, i3];
    roi_data = {Zslice, Xleft, Xright, Ytop, Ybottom, nVoxels};
    
    roi_template_str = fileread([ds.dirs.tpl filesep 'TBV' filesep 'ROIs_template2.roi']);
    fid = fopen(fn_result,'a');
    fprintf(fid,roi_template_str,roi_data{:});
    fprintf(fid,'%d\t%d\t%d\n', ind3D');
    fclose(fid);    
end
fprintf('Size xyz: %d %d %d\n',size(dv,1),size(dv,2),size(dv,3))
fprintf('Finished transforming .nii to .roi for TBV.\n >> Check coordinates within 96x96x52 <<\n')
save([outdir filesep ds.idstr '_tbv_voi' outdir(end) '.mat'],'xyz');












