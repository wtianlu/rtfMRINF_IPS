% WP3 rt-fMRI NF
% Lulu Wang
% Estimate GLM for each individual neurofeedback scan
% -------------------------------------------------------------------------
% EDITS:
% 2020.06.06 Use GLM to validate TBV; only need beta maps
% 2020.03.11 The big purge

%% ------------------------------------------------------------------------
function taskGLM(outdir,funcfiles,names,onsets,durations,regressors)

TR = 2;

spm('defaults', 'FMRI');
clear matlabbatch

matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(funcfiles);
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.name = names{1};
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = onsets;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = durations;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {regressors};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];%[1 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
spm_jobman('run',matlabbatch);

% matlabbatch{2}.spm.stats.review.spmmat = {[outdir '/SPM.mat']};
% matlabbatch{2}.spm.stats.review.display.matrix = 1;
% matlabbatch{2}.spm.stats.review.print = 'pdf';
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[outdir '/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run',matlabbatch);

% Create contrasts
clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat = {[outdir '/SPM.mat']};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'regulate+';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'regulate-';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;
spm_jobman('run',matlabbatch);

% Show results
clear matlabbatch
matlabbatch{1}.spm.stats.results.spmmat = {[outdir '/SPM.mat']};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{1}.spm.stats.results.conspec.extent = 40;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
% matlabbatch{1}.spm.stats.results.export{1}.pdf = true;
spm_jobman('run',matlabbatch);
