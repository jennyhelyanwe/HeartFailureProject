function [interpolated_pressure] = LinearInterpolatePressure(raw_pressure, d_p, d_m)
%% Interpolate pressure points at each image frame assuming linear interpolation between pressure data points.  
% Author: ZJW
% Date: 20 Feb 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interpolated_pressure = zeros(1, d_m);
interpolated_pressure(1) = raw_pressure(1);
interpolated_pressure(end) = raw_pressure(end);

    p_m_ratio = (d_p-1)/(d_m-1);
    for m = 2:d_m-1
        idx_interp = (m-1) * p_m_ratio + 1;
        idx_before = floor(idx_interp);
        idx_after = ceil(idx_interp);
        fraction = idx_interp - idx_before;
        % Weigh the contribution of the pressure data point before and
        % after the interpolation point using linear interpolation. 
        interpolated_pressure(m) = raw_pressure(idx_before)*(1-fraction) + raw_pressure(idx_after)*fraction;
    end

end