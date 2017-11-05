function [eIVR, ECG, a2a, n_c] = PWP(study_name, R)
%% Process pulmonary wedge pressure raw traces.  
% Read pulmonary wedge pressure trace and use this to identify the end of 
% isovolumic relaxation (and the beginning of passive filling). 
% Author: ZJW
% Date of first commit: 20 Feb 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in AOP trace and identify flagged points. 
files = dir('RawPlusECGTraces');
cd('RawPlusECGTraces');
files = files(~([files.isdir]));
num_files = length(files);
n = 1;
i_a = 1; % Counter for point of highest LA pressure. 
i_v = 1; % Counter for point of mitral valve opening.
i_r = 1; % Counter for R peaks of ECG.
i = 1; % Counter for pressure trace.
t_resolution = 0.004; % Temporal resolution of pressure trace is 4 ms.
while n <= num_files
    fname = files(n).name; 
    if ~isempty(strfind(fname, study_name))
        if strfind(fname, 'PW') % Read every PWP raw trace. 
            disp(fname); % Show current file being read to console. 
            fid = fopen(fname, 'r');
            header = fgetl(fid);
            idx = find(~cellfun(@isempty, strfind(strread(header, '%s', 'delimiter', '\t| '), 'PW'))); % Find column index for PW trace. 
            LV_column = idx(1); 
            ECG_columns = find(~cellfun(@isempty, strfind(strread(header, '%s', 'delimiter', '\t| '), 'I'))); % Find column index for I trace of ECG. 
            % Read PWP trace.             
            line = fgetl(fid);
            while ~feof(fid)
                values = strread(line, '%s', 'delimiter', '\t| ');
                temp = values{LV_column};
                % Convert value to string and check for flags
                if ~isempty(strfind(temp, 'A')) % 'A' flag marks one peak. 
                    A(i_a) = i;
                    At(i_a) = (i-1)*t_resolution;
                    i_a = i_a + 1;
                    l = length(temp);
                    p = temp(1:(l-3));
                    PWPressure(i) = str2num(p);
                elseif strfind(temp, 'V') % 'V' flag marks end IVR point. 
                    V(i_v) = i;
                    Vt(i_v) = (i-1)*t_resolution;
                    i_v = i_v + 1;
                    l = length(temp);
                    p = temp(1:(l-3));
                    PWPressure(i) = str2num(p);
                else
                    PWPressure(i) = str2num(temp);
                end
                t(i) = (i-1)*t_resolution; 
                % Get R peaks from ECG.
                temp_I = values{ECG_columns(1)};
                temp_II = values{ECG_columns(2)};
                if ~isempty(strfind(temp_I, 'R'))
                    R(i_r) = i;
                    Rt(i_r) = (i-1)*t_resolution;
                    i_r = i_r + 1;
                    l = length(temp_I);
                    v = temp_I(1:(l-3));
                    ECG_reading(i) = str2num(v);
                elseif ~isempty(strfind(temp_II, 'R'))
                    R(i_r) = i;
                    Rt(i_r) = (i-1)*t_resolution;
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
    n = n + 1; % Increment number of files processed. 
end
cd('../');      

%% Plot raw trace and landmark points
figure('units', 'normalized', 'outerposition', [0 0 1 1])
plot(t, PWPressure);
hold on
plot(At, PWPressure(A), 'r.', 'MarkerSize', 15); % Peak landmark.
plot(Vt, PWPressure(V), 'k.', 'MarkerSize', 15); % end IVR landmark. 
title('Raw Pulmonary Wedge Pressure trace');
ylabel('Pressure (mmHg)');   
xlabel('Time (s)');
for i = 1:length(R)
    plot(R(i)*t_resolution, 0:25, 'b.');
end

%% User select cycles using peak and eIVR landmarks. 
% Selected cycles are stored in HaemoAnalysis.mat. 
load('HaemoAnalysis');
[bool, idx] = ismember(study_name, HaemoAnalysis.study_name); 
if bool % If user selection has been done before, use the same cycles. 
    cycles = HaemoAnalysis.PWP{idx};
    n_c = length(cycles);
else % Get user input for cycle selection and save to .mat file for future use. 
    cycles = input('Enter array containing chosen cycle peak numbers (format [###]): ');
    n_c = length(cycles);
    if n_c == 0
        eIVR = [];
        ECG = [];
        a2a = [];
        n_c = 0;
        return
    end
    % Save cycles information in struct.
    HaemoAnalysis.PWP{end+1} = cycles;
    save('HaemoAnalysis', 'HaemoAnalysis'); % Save user selections to .mat file. 
end
% Initialise landmark indices arrays. 
eIVR.i = zeros(1, n_c);
ECG.i = zeros(1, n_c);
a2a.i = zeros(1, n_c);
a2a.t = zeros(1, n_c);
cycles_a = cycles;
cycles_v = zeros(1, n_c);
cycles_r = zeros(1, n_c);
% Find V waves and R peaks enclosed in selected cycle.
for i = 1:n_c
    for j = 1:length(V)
        if V(j) > A(cycles_a(i)) && V(j) < A(cycles_a(i)+1)
            cycles_v(i) = j;
            break
        end
    end
    for k = 1:length(R)
        if R(k) > A(cycles_a(i)) && R(k) < A(cycles_a(i)+1)
            cycles_r(i) = k;
            break
        end
    end
    a2a.i(i) = A(cycles_a(i)+1) -  A(cycles_a(i));
    a2a.t(i) = a2a.i(i)*t_resolution;
end
% Evaluate indices for end IVR and R peak with respect to where A peak
% occurs. 
eIVR.i = V(cycles_v) - A(cycles_a);
eIVR.t = eIVR.i*t_resolution;
ECG.i = R(cycles_r) - A(cycles_a);
ECG.t = ECG.i*t_resolution;
return