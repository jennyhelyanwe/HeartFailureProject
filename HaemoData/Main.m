%% Main script for running haemodynamic analysis for all studies. 
% Loop through all studies and perform haemodynamic analysis. 
% Author: ZJW
% Date of first commit: 20 Feb 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
study_names = {'STF_16', 'STF_17', 'STF_18','STF_19', 'STF_20',  'MR_250250'...
    , 'MR_293293', 'MR_119119', 'MR_054054', 'MR_104104', 'MR_236236', 'MR_269269','MR_087087', ...
    'MR_091091', 'MR_124124', 'MR_126126', 'STF_01', 'STF_02', 'STF_08', 'STF_09', 'STF_13', ...
    'MR_042042', 'STF_10', 'STF_11', 'MR_262262', 'STF_06', 'STF_12', 'MR_160160'};
for i = study_names
    % Run analysis for each study. 
    HaemoMain(char(i));
end

chdir('RegisteredPressure');
PassivePressureVolumeCurves;
EDP_Boxplots;
