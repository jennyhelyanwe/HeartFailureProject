function [R] = ECG(fname)
%% This function returns all the R peaks of the input ECG trace.

fid = fopen(fname);
I = fscanf(fid, '%s', 1);
II = fscanf(fid, '%s', 1);
i_a = 1; % Counter for point of highest LA pressure. 
i_v = 1; % Counter for point of mitral valve opening. 
i = 1; % Counter for pressure trace. 
while ~isempty(I)
    % Convert value to string and check for flags
    if ~isempty(strfind(I, 'R'))
        R(i_a) = i;
        Rt(i_a) = (i-1)*0.004;
        i_a = i_a + 1;
        l = length(I);
        p = I(1:(l-3));
        ECG(i) = str2num(p);
        t(i) = (i-1)*0.004;
    elseif ~isempty(strfind(II, 'R'))
        R(i_a) = i;
        Rt(i_a) = (i-1)*0.004;
        i_a = i_a + 1;
        l = length(II);
        p = II(1:(l-3));
        ECG(i) = str2num(p);
        t(i) = (i-1)*0.004;
    else
        ECG(i) = str2num(I);
        t(i) = (i-1)*0.004;
    end
    i = i + 1;
    I = fscanf(fid, '%s', 1);
    II = fscanf(fid, '%s', 1);
end


