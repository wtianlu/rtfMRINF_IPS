# rtfMRINF_IPS
 Modulating the interhemispheric activity balance in the intraparietal sulcus using real-time fMRI neurofeedback: development and proof-of-concept

## Overview 

This repository contains the scripts needed to reproduce the setup in the paper "Modulating the interhemispheric activity balance in the intraparietal sulcus using real-time fMRI neurofeedback: development and proof-of-concept".

## Contents

* Paradigm
  * Online NF code
    * main_nf.m: main experiment code
    * ds_check.m: script for online data quality checking
    * Templates/: NIFTI ROI and Turbo Brainvoyager settings templates to be customised for each individual participant 
    * functions/
  * CombiTVA assessment
    * CombiTVA.py: main experiment code
    * iViewXAPI.py: SMI eye tracker script
    * masks/
  * Questionnaires: TBA
* Analysis code
  * Functions/
  * WP3_main_analysis.m
  * WP3_1_prepfiles_nii.m
  * WP3_2_preproc_nii.m
  * WP3_3_analyse_nii.m
  * WP3_3_analyse_tva.m
  * WP3_4_analyse_fc.m
  * vis_questionnairedata.m
* Documentation: TBA

## Software requirements

The Matlab scripts have been tested in Matlab 2016b. The Python scripts have been tested in Python 2.7. 

## Citation

For usage of the scripts and the associated manuscript, please use the following:

Tianlu Wang, Ronald Peeters, Dante Mantini, Céline R. Gillebert (2020). Modulating the interhemispheric activity balance in the intraparietal sulcus using real-time fMRI neurofeedback: development and proof-of-concept. NeuroImage: Clinical. doi:10.1016/j.nicl.2020.102513
