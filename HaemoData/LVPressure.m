function [DS, ED, RR, LVP, t_all, n_c] = LVPressure(study_name)
%% Process left ventricular pressure raw traces. 
% Read in all LVP files for analysis and perform cycle selection. 
% Author: ZJW
% Date of first commit: 20 Feb 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in LVP trace and identify flagged points. 
files = dir('RawPlusECGTraces');
cd('RawPlusECGTraces');
files = files(~([files.isdir]));
num_files = length(files);
n = 1;
i_d = 1; % Counter for diastasis points.
i_s = 1; % Counter for peak pressure points.
i_e = 1; % Counter for end diastole points.
i_r = 1; % Counter for R peak of ECG. 
i = 1; % Counter for pressure trace.
t_resolution = 0.004; % Temporal resolution of pressure trace is 4 ms. 
while n <= num_files
    fname = files(n).name;
    if ~isempty(strfind(fname, study_name))
        if strfind(fname, 'LV') % Read every LV raw trace. 
            disp(fname); % Show current file being read to console. 
            fid = fopen(fname, 'r');
            header = fgetl(fid);
            idx = find(~cellfun(@isempty, strfind(strread(header, '%s', 'delimiter', '\t| '), 'LV')));
            LV_column = idx(1);
            ECG_columns = find(~cellfun(@isempty, strfind(strread(header, '%s', 'delimiter', '\t| '), 'I')));
            % Read LVP trace. 
            line = fgetl(fid);
            while ~feof(fid)
                values = strread(line, '%s', 'delimiter', '\t| ');
                temp = values{LV_column};
                % Convert value to string and check for flags
                if ~isempty(strfind(temp, 'S')) % 'S' flag marks peak pressure. 
                    S(i_s) = i;
                    St(i_s) = (i-1)*t_resolution;
                    i_s = i_s + 1;
                    l = length(temp);
                    p = temp(1:(l-3));
                    Pressure(i) = str2num(p);
                elseif strfind(temp, 'D') % 'D' flag marks lowest pressure (diastasis).
                    D(i_d) = i;
                    Dt(i_d) = (i-1)*t_resolution;
                    i_d = i_d + 1;
                    l = length(temp);
                    p = temp(1:(l-3));
                    Pressure(i) = str2num(p);
                elseif strfind(temp, 'E') % 'E' flag marks end of diastole. 
                    E(i_e) = i;
                    Et(i_e) = (i-1)*t_resolution;
                    i_e = i_e + 1;
                    l = length(temp);
                    p = temp(1:(l-3));
                    Pressure(i) = str2num(p);
                elseif strfind(temp, 'N/A')
                    Pressure(i) = 0;
                else
                    Pressure(i) = str2num(temp);
                end
                % Get R peaks. 
                temp_I = values{ECG_columns(1)};
                temp_II = values{ECG_columns(2)};
                if ~isempty(strfind(temp_I, 'R')) % 'R' peaks should coincide with ED. 
                    R(i_r) = i;
                    Rt(i_r) = (i-1)*0.004;
                    i_r = i_r + 1;
                    l = length(temp_I);
                    v = temp_I(1:(l-3));
                    ECG_reading(i) = str2num(v);
                elseif ~isempty(strfind(temp_II, 'R'))
                    R(i_r) = i;
                    Rt(i_r) = (i-1)*0.004;
                    i_r = i_r + 1;
                    l = length(temp_II);
                    v = temp_II(1:(l-3));
                    ECG_reading(i) = str2num(v);
                else
                    ECG_reading(i) = str2num(temp_I);
                end
                t(i) = (i-1)*t_resolution;
                line = fgetl(fid);
                i = i + 1;
            end
        end
    end
    n = n + 1;
end
cd('../');
%Filtering only for studies MR_119119 (N=45) and STF_11 (N=25): 
Pressure_prefilter = Pressure;
if strcmp(study_name, 'MR_119119')
    Pressure = medfilt1(Pressure, 45);
elseif strcmp(study_name, 'MR_152152')
    Pressure = medfilt1(Pressure, 25);
elseif strcmp(study_name, 'MR_042042')
    Pressure = medfilt1(Pressure, 25);
else
    Pressure = medfilt1(Pressure, 20); % Apply general filter to the trace to get rid of artefacts. 
end

%% Plot raw trace and landmark points. 
raw = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot(t, Pressure, 'k', 'LineWidth', 2);
hold on
plot(t, Pressure_prefilter, 'm', 'LineWidth', 0.7);
plot(St, Pressure(S), 'r.', 'MarkerSize', 30); % Plot peak points in red. 
plot(Dt, Pressure(D), 'b.', 'MarkerSize', 30); % Plot diastasis in blue. 
plot(Et, Pressure(E), 'g.', 'MarkerSize', 30); % Plot ED in green. 
plot(t(R), Pressure(R), 'bo', 'MarkerSize', 5); % Plot R peak occurence in blue circle. 
title('Pressure trace');
ylabel('LV pressure (mmHg)');
xlabel('Time (s)');
legend('Pressure trace', 'Pressure pre-filtering', 'Peaks', 'DS', 'ED', 'ECG R peak');
set(gca, 'fontsize', 16);
%axis([0 13 0 180]);
plot(t, repmat(15, 1, length(t)), 'r-');
print(raw, '-djpeg', ['AnalysisFigures/', study_name, '_LVP_raw']);

%% User select cycles using peak points
% Selected cycles are stored in HaemoAnalysis.mat. 
load('HaemoAnalysis');
[bool, idx] = ismember(study_name, HaemoAnalysis.study_name);
if bool % If user selection has been done before, use the same cycles. 
    cycles = HaemoAnalysis.LVP{idx};
    n_c = length(cycles);
else % Get user input for cycle selection and save to .mat file for future use. 
    cycles = input('Enter array containing chosen cycle peak numbers (format [###]): ');
    n_c = length(cycles);
    % Save cycles information in struct. 
    HaemoAnalysis.LVP{end+1} = cycles;
    save('HaemoAnalysis', 'HaemoAnalysis'); % Save user selections to .mat file. 
end
% Initialise landmark indices arrays. 
cycles_peak = cycles;
cycles_ed = zeros(1, n_c);
cycles_ds = zeros(1, n_c);
cycles_r = zeros(1, n_c);
% Find ED, diastasis, and peak points enclosed in selected cycle. 
for i = 1:n_c
    for j = 1:length(E)
        if E(j) > S(cycles_peak(i)) && E(j) < S(cycles_peak(i)+1)
            cycles_ed(i) = j;
        end
    end
    for k = 1:length(D)
        if D(k) > S(cycles_peak(i)) && D(k) < S(cycles_peak(i)+1)
            cycles_ds(i) = k;
        end
    end
    for l = 1:length(R)
        if R(l) > S(cycles_peak(i)) && R(l) < S(cycles_peak(i)+1)
            cycles_r(i) = l;
        end
    end
end
ECG.i = R(cycles_r);

%% For each cycle get sections from peak to DS, and from DS to ED. 
DS.i = zeros(1, n_c);
DS.t = zeros(1, n_c);
DS.p = zeros(1, n_c);
ES.i = zeros(1, n_c);
ES.t = zeros(1, n_c);
ES.p = zeros(1, n_c);
RR.i = zeros(1, n_c);
RR.t = zeros(1, n_c);
LVP = cell(n_c);
t_all = cell(n_c);
%ED_toggle = 0; % Use cath landmark to identify end diastole. 
ED_toggle = 1; % Use R peaks from ECG to identify end diastole. 
if ED_toggle == 1
    disp('Using R peaks from ECG to identify ED in each LVP cycle.');
else
    disp('Using cath. software calculated ED points in each LVP cycle.');
end
for j = 1:n_c
    % Get current entire cycle and convert to kPa.  
    LVP_c = Pressure(S(cycles_peak(j)):S(cycles_peak(j)+1))/7.50061561303;
    LVP{j} = LVP_c;
    % Get R to R interval of current cycle. 
    t_c = t(S(cycles_peak(j)):S(cycles_peak(j)+1)) - t(S(cycles_peak(j)));
    t_all{j} = t_c;
    RR.t(j) = t(S(cycles_peak(j)+1)) - t(S(cycles_peak(j)));
    RR.i(j) = S(cycles_peak(j)+1) - S(cycles_peak(j))+1;
    % Locate DS
    temp = find(LVP_c == min(LVP_c));
    DS.i(j) = temp(end);
    DS.t(j) = t_c(DS.i(j));
    DS.p(j) = LVP_c(DS.i(j));
    if ED_toggle == 1
        % Locate ED using R peak timing. 
        ED.i(j) = R(cycles_r(j)) - S(cycles_peak(j));
        ED.t(j) = t_c(ED.i(j));
        ED.p(j) = LVP_c(ED.i(j));
    else
        % Locate ED using machine indicated point. 
        ED.i(j) = E(cycles_ed(j)) - S(cycles_peak(j));
        ED.t(j) = t_c(ED.i(j));
        ED.p(j) = LVP_c(ED.i(j));
    end
end
return