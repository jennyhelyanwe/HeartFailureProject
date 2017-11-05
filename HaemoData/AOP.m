function [eIVC, ES, n_c] = AOP(study_name)
load('HaemoAnalysis');
%% Read in AOP trace and identify flagged points. 
% Read in raw traces. 
files = dir('RawPlusECGTraces');
cd('RawPlusECGTraces');
files = files(~([files.isdir]));
num_files = length(files);
n = 1;
i_s = 1;
i_e = 1;
i = 1;
disp('Reading in raw trace files...');
while n <= num_files
    fname = files(n).name;
    if ~isempty(strfind(fname, study_name))
        if strfind(fname, 'AO')
            disp(fname);
            fid = fopen(fname, 'r');
            header = fgetl(fid);
            idx = find(~cellfun(@isempty, strfind(strread(header, '%s', 'delimiter', '\t| '), 'AO')));
            column = idx(1);
            line = fgetl(fid);
            while ~feof(fid)
                values = strread(line, '%s', 'delimiter', '\t| ');
                temp = values{column};
                if strfind(temp, 'S')
                    S(i_s) = i;
                    St(i_s) = (i-1)*0.004; % Time stamp of ES.
                    i_s = i_s + 1;
                    l = length(temp);
                    p = temp(1:(l-3));  % Extract pressure value, get rid of string flag. 
                    Pressure(i) = str2num(p);
                elseif strfind(temp, 'D')
                    D(i_e) = i;
                    Dt(i_e) = (i-1)*0.004;
                    i_e = i_e + 1;
                    l = length(temp);
                    p = temp(1:(l-3));  % Extract pressure value, get rid fo string flag. 
                    Pressure(i) = str2num(p);
                elseif strfind(temp, 'N/A')
                    Pressure(i) = 0;
                else
                    Pressure(i) = str2num(temp);
                end
                t(i) = (i-1)*0.004;
                line = fgetl(fid);
                i = i + 1;
            end
        end
    end  
    n = n + 1;
end
cd('../')

%% Plot raw AOP trace with landmark points. 
figure('units', 'normalized', 'outerposition', [0 0 1 1])
plot(t, Pressure);
hold on
plot(St, Pressure(S), 'r.', 'MarkerSize', 15);
plot(Dt, Pressure(D), 'k.', 'MarkerSize', 15);
title('Raw aortic pressure trace');
ylabel('Pressure (mmHg)');
xlabel('Time (s)');

%% User select cycles using peak points. 
[bool, idx] = ismember(study_name, HaemoAnalysis.study_name);
if bool
    cycles = HaemoAnalysis.AOP{idx};
    n_c = length(cycles);
else
    cycles = input('Enter array containing chosen cycle peak numbers (format [###]): ');
    n_c = length(cycles);
    % Save cycles information in struct.
    HaemoAnalysis.AOP{end+1} = cycles; 
    save('HaemoAnalysis', 'HaemoAnalysis');
end

%% For each cycle, evaluate ES timing and end IVC timing as well as pressure
% values. 
eIVC.i = zeros(1, n_c);
eIVC.t = zeros(1, n_c);
eIVC.p = zeros(1, n_c);
ES.i = zeros(1, n_c);
ES.t = zeros(1, n_c);
ES.p = zeros(1, n_c);
for i = 1:n_c
    % Convert pressure to kPa.
    AOP_c = Pressure(S(cycles(i)):S(cycles(i)+1))*0.1333223899999367;
    RR_c = t(S(cycles(i)):S(cycles(i)+1)) - t(S(cycles(i))); % R-R interval of cycle. 
    % Locate end IVC.
    temp = find(AOP_c == min(AOP_c));
    eIVC.i(i) = temp(end);
    eIVC.t(i) = t(eIVC.i(i));
    eIVC.p(i) = AOP_c(eIVC.i(i));
    % Locate ES.
    AOP_c_dd = zeros(1, size(AOP_c, 2));
    for j = 2:(size(AOP_c, 2)-1)
        AOP_c_dd(j) = (AOP_c(j+1) - 2*AOP_c(j) + AOP_c(j-1))/(0.004)^2;
    end
    AOP_c_dd_section = AOP_c_dd(1:ceil(length(AOP_c_dd)/3));
    i_temp = find(AOP_c_dd_section==max(AOP_c_dd_section));
    ES.i(i) = i_temp(end);
    ES.t(i) = t(ES.i(i));
    ES.p(i) = AOP_c(ES.i(i));
end

return