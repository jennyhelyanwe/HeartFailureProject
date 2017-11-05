function [] = HaemoMain(study_name)
%% Haemodynamic analysis
% This script handles the analyses of the catheterisation data and
% quantifies the uncertainties related to both the timing and magnitude of
% pressure measurements.
% Author: ZJW
% Date of first commit: 20 Feb 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all;

%% Read in study name and details. 
% study_name = input('Enter study name: ', 's');
i = 0;
while i == 0
    i = 1;
    fid = fopen('StudyFrameNumbers.txt', 'r');
    temp = fgetl(fid);
    while isempty(strfind(temp, study_name))
        if temp == -1
            % Reached end of file without finding the study. 
            study_name = input('The study number you have entered is not found, please enter another: ', 's');
            i = 0;
            fclose(fid);
            break
        end
        temp = fgetl(fid);
    end
end
% Get important frame numbers and image paths from text file. 
temp = regexp(temp, '\s+', 'split');
study_name = temp{1};
mri.image_name = temp{2};
mri.num_frames = str2num(temp{3});
mri.ED = str2num(temp{4});
mri.eIVC = str2num(temp{5});
mri.ES = str2num(temp{6});
mri.eIVR = str2num(temp{7});
mri.DS = str2num(temp{8});
mri.volume_path = temp{9};

%% Process haemodynamic traces. 
haemo.t_resolution = 0.004; % Cath traces temporal resolution is 4 ms. 
% Analyse AOP
[haemo.AO_eIVC, haemo.AO_ES, haemo.n_AOP] = AOP(study_name);

% Analyse PWP. 
files = dir('RawPlusECGTraces');
num_files = length(files);
n = 1;
haemo.pwp_toggle = 0; % Default that PWP traces don't exist. 
while n <= num_files
    fname = files(n).name;
    if ~isempty(strfind(fname, study_name))
        if strfind(fname, 'PW')
            haemo.pwp_toggle = 1; % If found raw PWP trace turn on toggle. 
        end
    end
    n = n + 1;
end
if haemo.pwp_toggle
    % Analyse PWP
    [haemo.PW_eIVR, haemo.ECG_PWP, haemo.a2a, haemo.n_PWP] = PWP(study_name);
    if isempty(haemo.PW_eIVR)
        haemo.pwp_toggle = 0;
        load('HaemoAnalysis')
        [bool, ~] = ismember(study_name, HaemoAnalysis.study_name);
        if bool == 0
            % Save no cycles for this study in the struct. 
            HaemoAnalysis.PWP{end+1} = [];
            save('HaemoAnalysis', 'HaemoAnalysis')
        end
    else
        haemo.pwp_toggle = 1;
    end
else
    load('HaemoAnalysis')
    [bool, ~] = ismember(study_name, HaemoAnalysis.study_name);
    if bool == 0
        % Save no cycles for this study in the struct. 
        HaemoAnalysis.PWP{end+1} = [];
        save('HaemoAnalysis', 'HaemoAnalysis')
    end
    haemo.pwp_toggle = 0;
end

% Analyse LVP
[haemo.DS, haemo.ED, haemo.RR, haemo.LVP_cycles, haemo.t_cycles, haemo.n_LVP] = LVPressure(study_name);

% Save cycles information in structure. 
load('HaemoAnalysis')
[bool, ~] = ismember(study_name, HaemoAnalysis.study_name);
if bool == 0
    HaemoAnalysis.study_name{end+1} = study_name;
    save('HaemoAnalysis', 'HaemoAnalysis');
end

%% Process cardiac MRI outputs. 
% Get image landmark timing. 
if exist(['ImageTiming/', study_name, '_MRI_info.mat'])
    load(['ImageTiming/', study_name, '_MRI_info.mat']);
    r_r = MRI_Info.RR_mean/1000; % R to R interval from each image. 
    int = MRI_Info.TS_mean/1000; % Time interval between each image. 
else
    cd('../HaemoData_NoUncertaintyQuantification/STF_Images/');
    [r_r] = ExtractTriggerTime(image_name, mri.num_frames, n);
    int = r_r/(mri.num_frames-1)/1000;
end
fprintf('The mean R-R interval of MRI is %.2fs. The frame-to-frame interval is %.2fs\n', r_r, int);
% Timing data for the MR images. 
mri.t = 0:int:r_r;
% Get LV volumes from CIM analysis. 
fname = ['../CIM_Models/Studies/', mri.volume_path];
fid = fopen(fname, 'r');
junk = fgets(fid);
i = 1;
while isempty(strfind(junk, 'TOTAL'))
    junk = fgets(fid);
    i = i + 1;
end
for i = 1:2
    junk = fgets(fid);
end
v_name = fscanf(fid, '%s', 2);
V = fscanf(fid, '%f', mri.num_frames);
mri.V = V';
fclose(fid);

%% Interpolate pressure at each MRI frame. 
disp('Aligning pressure with MR images...');
[output] = AlignPressureMRI(haemo, mri);

%% Write pressure values to text files.
disp('Writing pressure values to text files...');
% Write offset LVP to text file. 
filename = ['RegisteredPressure/', output.simulation_name, '_registered_LVP.txt'];
WriteToTextFile(output.LVP, filename, haemo.pwp_toggle); 
% Write no offset LV pressures to text file.
filename = ['RegisteredPressure/', output.simulation_name, '_registered_LVP_no_offset.txt'];
WriteToTextFile(output.LVP_no_offset, filename, haemo.pwp_toggle);
% Write offset pressure average
fw_mean = fopen(['RegisteredPressure/', output.simulation_name, '_registered_LVP_mean.txt'], 'w');
for i = 1:length(output.LVP_average)
    fprintf(fw_mean, '%d\t%f\t%f\n', i, output.LVP_average(i), output.LVP_ste(i));
end
% Write no offset pressure average
fw_no_mean = fopen(['RegisteredPressure/', output.simulation_name, '_registered_LVP_no_offset_mean.txt'], 'w');
for i = 1:length(output.LVP_average)
    fprintf(fw_no_mean, '%d\t%f\t%f\n', i, output.LVP_no_offset_average(i), output.LVP_no_offset_ste(i));
end
fclose(fw_mean);
fclose(fw_no_mean);

%% Generate extra figures for paper
% Plot passive pressure trace only. 
% figure
% hold on 
% for i = 1:n_LVP
%     plot(0:0.023:8*0.23, LVP.ds_ed(:, i), 'b*');
% end
% plot(0:0.023:8*0.23, [LVP_average(23:30), LVP_average(1)], 'k-');
% title('Passive LV pressure')
% xlabel('Time (s)')
% ylabel('LV pressure (kPa)')
% set(gca, 'fontsize', 16);

%% Save analysis outputs as .mat files. 
save(['RegisteredPressure/', output.simulation_name, '_mri'], 'mri');
save(['RegisteredPressure/', output.simulation_name, '_lvp'], 'output');

%% Wrap up. 
disp('Finished haemodynamic analysis.');