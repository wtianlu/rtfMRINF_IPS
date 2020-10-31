% 190730 integrate with ds
% 190626 nii2tbvroi new version

function nii2tbvroi(ds)
fprintf('***nii2tbvroi.m****************************\n')

% Copy template 1 
fn_result = [ds.dirs.rois filesep ds.idstr '_tbv_voi.roi'];
if exist(fn_result,'file')==0
    copyfile([ds.dirs.tpl filesep 'TBV' filesep 'ROIs_template1.roi'],fn_result);
else
    fprintf('%s already exists!\n',fn_result)
    return;
end

fp_roi = cellfun(@(x) [ds.dirs.rois filesep 'brwvoi' num2str(x) '.nii'],num2cell(1:2),'UniformOutput',false)';
for ii = 1:length(fp_roi)
    disp(fp_roi{ii})
    % load nifti
    V = spm_vol(fp_roi{ii});
    [dv,~] = spm_read_vols(V);
    % get coordinates
    ind=find(dv>0.8);
    [i1, i2, i3] = ind2sub(size(dv), ind); xyz{ii} = [i1, i2, i3]';
    nVoxels = length(i1);
    fprintf(' nVoxels: %d x: %d-%d y: %d-%d slices: %d-%d\n',nVoxels,min(i1),max(i1),min(i2),max(i2),min(i3),max(i3))
    
    % Get header info
    Zslice = round(mean(i3));
    Xleft = min(i1(i3==Zslice));
    Xright = max(i1(i3==Zslice)); 
    Ytop = min(i2(i3==Zslice));
    Ybottom = max(i2(i3==Zslice)); 
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
save([ds.dirs.rois filesep ds.idstr '_tbv_voi.mat'],'xyz');
fprintf('*******************************************\n')

