% PROJECT:      WP3 - rt-fMRI NF for self-regulation of interhemispheric IPS activity
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Analysis of all NIFTI data across all participants
% -------------------------------------------------------------------------
% 1. Task-based: estimate GLM per run
% -------------------------------------------------------------------------
% 2020.06.05 Stay in Native space, folder: run-1n, file: s4D.nii - ver1
% 2020.06.06 Fix co-registration. File: sr4D.nii - ver1_1

function dirs = WP3_3_analyse_nii(dirs,qc)
if nargin<2, qc=false; end

%% Initialisation
np = dirs.n.p; ns = dirs.n.s; nr = dirs.n.r; vols = 210;
curcodedir = '/Users/Lulu/Documents/Experiments/WP3a_rtfMRI_NF/Code/Final_200625';
%             cd('/Volumes/LACIE SHARE/WP3/Code')
cd(curcodedir)
load('visualisationsettings.mat')

fprintf('\n\n***************** Start WP3_3_analyse_nii *****************\n\n')
%% 1. GLM analysis - Run by run

for p = 1%:np
    for s = 1:ns
        for r = 1:nr
            cd(curcodedir)
            % prepare functional data and directory
            dir_out = sprintf('%ssub-0%d/ses-%d/run-%dn/',dirs.results.main,p,s,r); % edit 200625 2 regressors
            if exist(dir_out,'dir')==0, mkdir(dir_out); end
            
            dirs.results.func{r,s,p} = dir_out;

            funcfiles = cellfun(@(x) sprintf('%s,%d',strrep(dirs.processed.func{r,s,p},'sw','sr'),x),...
                num2cell(1:vols)','UniformOutput',false);

            % prepare conditions
            fn_cond = sprintf('%s/conditions.mat',dirs.data.func{r,s,p});
            load(fn_cond)%,'names','onsets','durations')

            % Check ART 
            load(dirs.data.artreg{r,s,p},'R')
            R = R(:,1:end-1); % remove framewise displacement

            % prepare regressors
            fn_reg = sprintf('%sreg_%s.txt',dirs.data.func{r,s,p},dirs.data.funcfn(1:end-4));
            if exist(fn_reg,'file')==0
                fid = fopen(fn_reg,'w');
                fprintf(fid,[repmat('%.6e ',1,size(R,2)-1) '%.6e\n'],R');
                fclose(fid);
            end
            dirs.results.regressors{r,s,p} = fn_reg;
            fprintf('prepared %s\n',fn_reg)

            if exist([dir_out 'SPM.mat'],'file')==0
                % estimate GLM
                taskGLM(dir_out,funcfiles,names,onsets,durations,fn_reg)
            fprintf('%sSPM.mat finished!\n',dir_out)
            else
                
            fprintf('%sSPM.mat already fitted!\n',dir_out)
            end
      
        end
    end
end


fprintf('\n\n***************** End WP3_3_analyse_nii *****************\n\n')

% errorspmfn = dir('/Volumes/LACIE SHARE/WP3/DataAnalysis/sub-0*/ses-*/run-*n/SPM.mat')
% for i = 1:length(errorspmfn)
%     delete([errorspmfn(i).folder '/' errorspmfn(i).name])
% end

