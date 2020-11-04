% PROJECT:      WP3 - rt-fMRI NF for self-regulation of interhemispheric IPS activity
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Main organisation script for preprocessing, analysis, visualisation
%               1. Check raw data
%               2. Pre-processing
%               3. Statistical analyses
%               4. Visualisation
% -------------------------------------------------------------------------
% EDITS:
% 2020.06.25 Skipping first 10 volumes
% 2020.06.05 Version without beta estimates and FC comparisons
% 2020.05.28 Clean and final version
% 2020.03.11 Initialisation
% -------------------------------------------------------------------------
%% Initialisation
clc; clearvars; close all;
curcodedir = '/Users/Lulu/Documents/Experiments/WP3a_rtfMRI_NF/Code/Final_200625/';

dirsmat = [curcodedir 'dirsmat.mat'];
if exist(dirsmat,'file')==0
%     dirs.main = '/Volumes/LACIE SHARE/WP3/';
    dirs.main = '/Users/Lulu/Documents/Experiments/WP3a_rtfMRI_NF/';
    dirs.raw.main = [dirs.main 'DataRaw/'];
    dirs.data.main = [dirs.main 'DataDerived/'];
    dirs.results.main = [dirs.main 'DataAnalysis/'];
    dirs.data.quest = [dirs.data.main 'Group/WP3_participants.txt'];
    dirs.n.p = 6; dirs.n.s = 3; dirs.n.r = 5; dirs.n.rest = 2;
else
    load(dirsmat)
end
fprintf('Initialisation in main folder: %s\n',dirs.main)
cd(fileparts(dirsmat))
qc = false; % quality check - set to true to check nii with SPM
fprintf('Ready!\n')
%% If finished, save dirsmat
save(dirsmat,'dirs')

%% 1. Prepare NIFTI files -------------------------------------------------
clc;
% Check NIFTI data
dirs = WP3_1_prepfiles_nii(dirs);
save(dirsmat,'dirs')
% Preprocess NIFTI data
dirs = WP3_2_preproc_nii(dirs);
cd(curcodedir)
save(dirsmat,'dirs')
% Fit NIFTI data
dirs = WP3_3_analyse_nii(dirs);

%% 2. Read betas + analysis + visualisation
dirs = WP3_4_analyse_betas(dirs);

% Read online dPSC data + analysis + visualisation
% dirs = WP3_4_analyse_dPSC(dirs);
dirs = WP3_4_analyse_dPSCbetas(dirs);

% Fit FC data + analysis + visualisation
dirs = WP3_3_analyse_fc(dirs);

% Analyse behavioural data
dirs = WP3_3_analyse_tva(dirs,qc);

% Visualise and analyse questionnaire info
vis_questionnairedata(dirs);





