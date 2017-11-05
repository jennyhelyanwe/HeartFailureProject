function [output] = AlignPressureMRI(haemo, mri)
%% Interpolate LV pressure at each MRI frame for current study.
% Author: ZJW
% Date: 20 Feb 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot all LVP raw traces in sections according to landmarks.
tracesFig = figure;
hold on
DSt_mean = mean(haemo.DS.t);
EDt_mean = mean(haemo.ED.t);
RRt_mean = mean(haemo.RR.t);
LVP_times = cell(haemo.n_LVP);
for i = 1:haemo.n_LVP
    t_peak2ds_scaled = linspace(0, DSt_mean, length(1:haemo.DS.i(i)));
    t_ds2ed_scaled = linspace(DSt_mean, EDt_mean, length(haemo.DS.i(i):haemo.ED.i(i)));
    t_ed2peak_scaled = linspace(EDt_mean, RRt_mean , length(haemo.ED.i(i):haemo.RR.i(i)));
    LVP_times{i} = [t_peak2ds_scaled, t_ds2ed_scaled(2:end), t_ed2peak_scaled(2:end)];
    figure(tracesFig);
    plot(LVP_times{i}, haemo.LVP_cycles{i}, 'LineWidth', 2);
end
plot([DSt_mean, DSt_mean], [0 25], 'b-', 'LineWidth', 2);
plot([EDt_mean, EDt_mean], [0 25], 'g-', 'LineWidth', 2);
title('Temporally scaled cardiac cycles');
set(gca, 'fontsize', 16);
axis([0 0.8 0 22])
xlabel('Time (s)')
ylabel('LV pressure (kPa)')
box on

%% Initialise figures to show analysis results.
imageFig = figure; % This will show the offset LVP traces.
hold on
imageNoOffsetFig = figure; % This will show the raw LVP traces with no offset.
hold on
PVFig = figure; % This will show the pressure volume curves.
hold on

%% Interpolate between DS and ED, this should have n_LVP alternatives.
disp('Interpolating between DS and ED...');
d_m_ds2ed = mri.num_frames - mri.DS + mri.ED + 1; % Number of frames in images from DS to ED.
d_p_ds2ed = haemo.ED.i - haemo.DS.i; % Number of data points in pressure from DS to ED.
LVP.ds_ed = zeros(d_m_ds2ed, haemo.n_LVP);
LV_no_offset.ds_ed = zeros(d_m_ds2ed, haemo.n_LVP);
for i = 1:haemo.n_LVP
    % Get DS to ED portion of current pressure cycle.
    LVP_c = haemo.LVP_cycles{i};
    ds_ed_p_raw = LVP_c(haemo.DS.i(i):haemo.ED.i(i)); % Pressure from ds to ed.
    % Interpolate pressure segment for each MRI frame.
    interpolated_pressure = LinearInterpolatePressure(ds_ed_p_raw, d_p_ds2ed(i), d_m_ds2ed);
    LVP_no_offset.ds_ed(:,i) = interpolated_pressure;
    LVP.ds_ed(:, i) = interpolated_pressure - haemo.DS.p(i);
    figure(imageFig);
    plot([mri.t(mri.DS:mri.num_frames), mri.t(mri.ED)], LVP.ds_ed(:, i), 'b*');
    figure(imageNoOffsetFig);
    plot([mri.t(mri.DS:mri.num_frames), mri.t(mri.ED)], LVP_no_offset.ds_ed(:, i), 'b*');
    figure(PVFig);
    plot([mri.V(mri.DS:mri.num_frames), mri.V(mri.ED)], LVP.ds_ed(:, i), 'b*');
end

%% Interpolation between ED & end IVC and between endIVC & ES has n_LVP*n_AOP alternatives.
disp('Interpolating between ED and end of IVC and then to ES...');
% Aortic inflection point marks end IVC, and splits the LVP trace from ED to ES
% into two portions.
figure(tracesFig);
for i = 1:haemo.n_AOP
    plot([0 1.5], [haemo.AO_ES.p(i), haemo.AO_ES.p(i)], 'r-');
    plot([0 1.5], [haemo.AO_eIVC.p(i), haemo.AO_eIVC.p(i)], 'm-');
end
d_m_ed2eivc = mri.eIVC - mri.ED + 1; % Number of frames in images from ED to end IVC.
d_m_eivc2es = mri.ES - mri.eIVC + 1; % Number of frames in images from end IVC to ES.
% Interpoalte pressure points at each image frame linearly.
LVP.ed_ievc = zeros(d_m_ed2eivc-1, haemo.n_LVP*haemo.n_AOP);
LVP_no_offset.ed_ievc = zeros(d_m_ed2eivc-1, haemo.n_LVP*haemo.n_AOP);
LVP.eivc_es = zeros(d_m_eivc2es-1, haemo.n_LVP*haemo.n_AOP);
LVP_no_offset.eivc_es = zeros(d_m_eivc2es-1, haemo.n_LVP*haemo.n_AOP);
count = 1;
eIVC.i = zeros(1, haemo.n_LVP*haemo.n_AOP); % Indices for end IVC landmarks.
ES.i = zeros(1, haemo.n_LVP*haemo.n_AOP); % Indices for ES landmarks.
for i = 1:haemo.n_LVP
    LVP_c = haemo.LVP_cycles{i}; % Get current LVP cycle.
    LVP_ed_peak = LVP_c(haemo.ED.i(i):end); % Cut out section of LVP from ED to peak.
    LVP_peak_ds = LVP_c(1:haemo.DS.i(i)); % Cut out section of LVP from peak to diastasis.
    for j = 1:haemo.n_AOP
        % Find end IVC using aortic pressure.
        temp = find(abs(LVP_ed_peak-haemo.AO_eIVC.p(j)) == min(abs(LVP_ed_peak - haemo.AO_eIVC.p(j)))); % Use AOP to find end IVC.
        eIVC.i(count) = temp(end)+haemo.ED.i(i);
        d_p_ed_eivc = eIVC.i(count) - haemo.ED.i(i); % Number of pressure data points from ED to end IVC.
        ed_eivc_raw = LVP_c(haemo.ED.i(i):eIVC.i(count)); % Pressure from ED to eIVC.
        % Interpolate pressure segment at each MRI frame.
        interpolated_pressure = LinearInterpolatePressure(ed_eivc_raw, d_p_ed_eivc, d_m_ed2eivc);
        LVP_no_offset.ed_eivc(:, count) = interpolated_pressure(2:end);
        LVP.ed_eivc(:, count) = interpolated_pressure(2:end) - haemo.DS.p(i);
        % Plot for figures.
        figure(imageFig);
        plot(mri.t((mri.ED+1):mri.eIVC), LVP.ed_eivc(:, count), 'g*');
        figure(imageNoOffsetFig);
        plot(mri.t((mri.ED+1):mri.eIVC), LVP_no_offset.ed_eivc(:, count), 'g*');
        figure(PVFig);
        plot(mri.V((mri.ED+1):mri.eIVC), LVP.ed_eivc(:, count), 'g*');
        
        % Interpolate pressure at images between end IVC to ES.
        temp = find(abs(LVP_peak_ds-haemo.AO_ES.p(j))==min(abs(LVP_peak_ds-haemo.AO_ES.p(j)))); % Use AOP to find ES.
        ES.i(count) = temp(end);
        ES.t(count) = ES.i(count)*haemo.t_resolution*haemo.DS.t(i)/mean(haemo.DS.t);
        d_p_eivc2es = ES.i(count) + length(LVP_c) - eIVC.i(count); % Number of pressure data points between end IVC and ES.
        eivc_es_raw = LVP_c([eIVC.i(count):end, 1:ES.i(count)]);
        % Interpolate pressure segment at each MRI frame.
        interpolated_pressure = LinearInterpolatePressure(eivc_es_raw, d_p_eivc2es, d_m_eivc2es);
        LVP_no_offset.eivc_es(:, count) = interpolated_pressure(2:end);
        LVP.eivc_es(:, count) = interpolated_pressure(2:end) - haemo.DS.p(i);
        % Plot analysis.
        figure(imageFig);
        plot(mri.t((mri.eIVC+1):mri.ES), LVP.eivc_es(:, count), 'r*');
        figure(imageNoOffsetFig);
        plot(mri.t((mri.eIVC+1):mri.ES), LVP_no_offset.eivc_es(:, count), 'r*');
        figure(PVFig);
        plot(mri.V((mri.eIVC+1):mri.ES), LVP.eivc_es(:, count), 'r*');
        count = count + 1;
    end
end

%% Interpolation between ES & end IVR and between end IVR & DS has n_LVP*n_AOP*n_PWP alternatives.
disp('Interpolating between ES and end of IVR and then to DS...');
if haemo.pwp_toggle == 1
    % PWP available, use V waves to locate end IVR.
    d_m_es2eivr = mri.eIVR - mri.ES + 1;
    LVP.es_eivr = zeros(d_m_es2eivr-1, haemo.n_LVP*haemo.n_AOP*haemo.n_PWP);
    LVP_no_offset.es_eivr = zeros(d_m_es2eivr-1, haemo.n_LVP*haemo.n_AOP*haemo.n_PWP);
    d_m_eivr2ds = mri.DS - mri.eIVR + 1;
    LVP.eivr_ds = zeros(d_m_eivr2ds-2, haemo.n_LVP*haemo.n_AOP*haemo.n_PWP);
    LVP_no_offset.eivr_ds = zeros(d_m_eivr2ds-2, haemo.n_LVP*haemo.n_AOP*haemo.n_PWP);
    es_count = 1;
    total_count = 1;
    used_pwp_count = 1;
    eIVR.i = zeros(1, haemo.n_LVP*haemo.n_AOP*haemo.n_PWP);
    eIVR.t = zeros(1, haemo.n_LVP*haemo.n_AOP*haemo.n_PWP);
    for i = 1:haemo.n_LVP
        LVP_c = haemo.LVP_cycles{i};
        for j = 1:haemo.n_AOP
            for k = 1:haemo.n_PWP
                % Find end IVR from PWP trace.
                PWP_scale_t = haemo.a2a.t(k)/haemo.RR.t(i);
                eIVR_temp = EDt_mean-(haemo.ECG_PWP.t(k)-haemo.PW_eIVR.t(k))*haemo.RR.t(i)/haemo.a2a.t(k);
                eIVR_temp = eIVR_temp*mean(haemo.DS.t)/haemo.DS.t(i);
                if (haemo.DS.t(i) - eIVR_temp <= 0.05*haemo.RR.t(i)) || (eIVR_temp - ES.t(es_count) <= 0.05*haemo.RR.t(i))
                    total_count = total_count + 1;
                    continue
                end
                eIVR.t(total_count) = eIVR_temp;
                eIVR.i(total_count) = find(abs(LVP_times{i}-eIVR.t(total_count))==min(abs(LVP_times{i}-eIVR.t(total_count))));
                % Interpolate pressure from ES to end IVR.
                d_p_es2eivr = eIVR.i(total_count) - ES.i(es_count);
                es_eivr_raw = LVP_c(ES.i(es_count):eIVR.i(total_count));
                interpolated_pressure = LinearInterpolatePressure(es_eivr_raw, d_p_es2eivr, d_m_es2eivr);
                LVP_no_offset.es_eivr(:, total_count) = interpolated_pressure(2:end);
                LVP.es_eivr(:, total_count) = interpolated_pressure(2:end) - haemo.DS.p(i);
                % Plot analysis.
                figure(imageFig);
                plot(mri.t((mri.ES+1):mri.eIVR), LVP.es_eivr(:,total_count), 'm*');
                figure(imageNoOffsetFig);
                plot(mri.t((mri.ES+1):mri.eIVR), LVP_no_offset.es_eivr(:,total_count), 'm*');
                figure(PVFig);
                plot(mri.V((mri.ES+1):mri.eIVR), LVP.es_eivr(:,total_count), 'm*');
                
                % end IVR to DS
                d_p_eivr2ds = haemo.DS.i(i) - eIVR.i(total_count);
                %if d_p_eivr2ds > d_m_eivr2ds
                eivr_ds_raw = LVP_c(eIVR.i(total_count):haemo.DS.i(i));
                interpolated_pressure = LinearInterpolatePressure(eivr_ds_raw, d_p_eivr2ds, d_m_eivr2ds);
                LVP_no_offset.eivr_ds(:, total_count) = interpolated_pressure(2:end-1);
                LVP.eivr_ds(:, total_count) = interpolated_pressure(2:end-1) - haemo.DS.p(i);
                % Plot analysis.
                figure(imageFig);
                plot(mri.t((mri.eIVR+1):(mri.DS-1)), LVP.eivr_ds(:, total_count), 'c*');
                figure(imageNoOffsetFig);
                plot(mri.t((mri.eIVR+1):(mri.DS-1)), LVP_no_offset.eivr_ds(:, total_count), 'c*');
                figure(PVFig);
                plot(mri.V((mri.eIVR+1):(mri.DS-1)), LVP.eivr_ds(:, total_count), 'c*');
                
                total_count = total_count + 1;
                %end
            end
            es_count = es_count + 1;
        end
    end
    eIVR.t(:, ~any(eIVR.t, 1)) = [];
    eIVR.i(:, ~any(eIVR.i, 1)) = [];
    LVP.es_eivr(:, ~any(LVP.es_eivr,1)) = [];
    LVP.eivr_ds(:, ~any(LVP.eivr_ds,1)) = [];
    LVP_no_offset.es_eivr(:, ~any(LVP_no_offset.es_eivr,1)) = [];
    LVP_no_offset.eivr_ds(:, ~any(LVP_no_offset.eivr_ds,1)) = [];
else
    % PWP is not available, then just linearly interpolate between ES to DS.
    d_m_es2ds = mri.DS - mri.ES + 1;
    LVP.es_ds = zeros(d_m_es2ds-2, haemo.n_LVP*haemo.n_AOP);
    LVP_no_offset.es_ds = zeros(d_m_es2ds-2, haemo.n_LVP*haemo.n_AOP);
    es_count = 1;
    for i = 1:haemo.n_LVP
        LVP_c = haemo.LVP_cycles{i};
        t_c = haemo.t_cycles{i};
        for j = 1:haemo.n_AOP
            d_p_es2ds = haemo.DS.i(i) - ES.i(es_count);
            es_ds_raw = LVP_c(ES.i(es_count):haemo.DS.i(i));
            interpolated_pressure = LinearInterpolatePressure(es_ds_raw, d_p_es2ds, d_m_es2ds);
            LVP_no_offset.es_ds(:, es_count) = interpolated_pressure(2:end-1);
            LVP.es_ds(:, es_count) = interpolated_pressure(2:end-1) - haemo.DS.p(i);
            % Plot analysis.
            figure(imageFig);
            plot(mri.t((mri.ES+1):(mri.DS-1)), LVP.es_ds(:, es_count), 'm*');
            figure(imageNoOffsetFig);
            plot(mri.t((mri.ES+1):(mri.DS-1)), LVP_no_offset.es_ds(:,es_count), 'm*');
            figure(PVFig);
            plot(mri.V((mri.ES+1):(mri.DS-1)), LVP.es_ds(:,es_count), 'm*');
            es_count = es_count + 1;
        end
    end
    %     %so use 2nd derivative of LVP trace to locate
    %     % the end of IVR.
    %     d_m_es2eivr = mri.eIVR - mri.ES + 1;
    %     LVP.es_eivr = zeros(d_m_es2eivr-1, haemo.n_LVP*haemo.n_AOP);
    %     LVP_no_offset.es_eivr = zeros(d_m_es2eivr-1, haemo.n_LVP*haemo.n_AOP);
    %     d_m_eivr2ds = mri.DS - mri.eIVR + 1 ;
    %     LVP.eivr_ds = zeros(d_m_eivr2ds-2, haemo.n_LVP*haemo.n_AOP);
    %     LVP_no_offset.eivr_ds = zeros(d_m_eivr2ds-2, haemo.n_LVP*haemo.n_AOP);
    %     eIVR.i = zeros(1, haemo.n_LVP*haemo.n_AOP);
    %     eIVR.t = eIVR.i;
    %     eIVR_count = 1;
    %     for i = 1:haemo.n_LVP
    %         LVP_c = haemo.LVP_cycles{i};
    %         t_c = haemo.t_cycles{i};
    %         for j = 1:haemo.n_AOP
    %             LVP_es_ds = LVP_c(ES.i(eIVR_count):haemo.DS.i(i));
    %             % Use 2nd derivative of LVP to identify end IVR.
    %             LVP_dd = zeros(1, length(LVP_es_ds));
    %             for k = 2:(size(LVP_es_ds, 2)-1)
    %                 LVP_dd(k) = (LVP_es_ds(k+1) - 2*LVP_es_ds(k) + LVP_es_ds(k-1))/(haemo.t_resolution)^2;
    %             end
    %             [~, locs] = findpeaks(LVP_dd, 'SortStr', 'descend');
    %             if ~isempty(locs)
    %                 eIVR.i(eIVR_count) = locs(end) + ES.i(eIVR_count);
    %             end
    %             d_p_eivr2ds = haemo.DS.i(i) - eIVR.i(eIVR_count);
    %             d_p_es2eivr = eIVR.i(eIVR_count) - ES.i(eIVR_count);
    %             if (d_p_eivr2ds > d_m_eivr2ds) && (d_p_es2eivr > d_m_es2eivr)
    %                 eIVR.t(eIVR_count) = t_c(eIVR.i(eIVR_count));
    %                 % Interpolate pressure from ES to end IVR.
    %                 es_eivr_raw = LVP_c(ES.i(eIVR_count):eIVR.i(eIVR_count));
    %                 interpolated_pressure = LinearInterpolatePressure(es_eivr_raw, d_p_es2eivr, d_m_es2eivr);
    %                 LVP_no_offset.es_eivr(:, eIVR_count) = interpolated_pressure(2:end);
    %                 LVP.es_eivr(:, eIVR_count) = interpolated_pressure(2:end) - haemo.DS.p(i);
    %                 % Plot analysis.
    %                 figure(imageFig);
    %                 plot(mri.t((mri.ES+1):mri.eIVR), LVP.es_eivr(:,eIVR_count), 'm*');
    %                 figure(imageNoOffsetFig);
    %                 plot(mri.t((mri.ES+1):mri.eIVR), LVP_no_offset.es_eivr(:,eIVR_count), 'm*');
    %                 figure(PVFig);
    %                 plot(mri.V((mri.ES+1):mri.eIVR), LVP.es_eivr(:,eIVR_count), 'm*');
    %
    %                 % end IVR to DS
    %                 eivr_ds_raw = LVP_c(eIVR.i(eIVR_count):haemo.DS.i(i));
    %                 interpolated_pressure = LinearInterpolatePressure(eivr_ds_raw, d_p_eivr2ds, d_m_eivr2ds);
    %                 LVP_no_offset.eivr_ds(:, eIVR_count) = interpolated_pressure(2:end-1);
    %                 LVP.eivr_ds(:, eIVR_count) = interpolated_pressure(2:end-1) - haemo.DS.p(i);
    %                 % Plot analysis.
    %                 figure(imageFig);
    %                 plot(mri.t((mri.eIVR+1):(mri.DS-1)), LVP.eivr_ds(:, eIVR_count), 'c*');
    %                 figure(imageNoOffsetFig);
    %                 plot(mri.t((mri.eIVR+1):(mri.DS-1)), LVP_no_offset.eivr_ds(:, eIVR_count), 'c*');
    %                 figure(PVFig);
    %                 plot(mri.V((mri.eIVR+1):(mri.DS-1)), LVP.eivr_ds(:, eIVR_count), 'c*');
    %                 eIVR_count = eIVR_count + 1;
    %             end
    %         end
    %     end
end
%% Put all interpolated segments of pressure together.
disp('Concatenating pressure segments...');
output.LVP = LVP;
output.LVP_no_offset = LVP_no_offset;
if haemo.pwp_toggle == 1
    output.LVP_average = [mean(LVP.ds_ed(end,:)); mean(LVP.ed_eivc, 2);
        mean(LVP.eivc_es, 2); mean(LVP.es_eivr, 2); mean(LVP.eivr_ds, 2);
        mean(LVP.ds_ed(1:end-1, :), 2)];
    output.LVP_ste = [std(LVP.ds_ed(end,:))/sqrt(length(LVP.ds_ed(end,:)));
        std(LVP.ed_eivc,0,2)/sqrt(length(LVP.ed_eivc));
        std(LVP.eivc_es,0,2)/sqrt(length(LVP.eivc_es));
        std(LVP.es_eivr,0,2)/sqrt(length(LVP.es_eivr));
        std(LVP.eivr_ds,0,2)/sqrt(length(LVP.eivr_ds));
        std(LVP.ds_ed(1:end-1, :),0,2)/sqrt(length(LVP.ds_ed(1:end-1, :)))];
    output.LVP_no_offset_average = [mean(LVP_no_offset.ds_ed(end,:));
        mean(LVP_no_offset.ed_eivc, 2); mean(LVP_no_offset.eivc_es, 2);
        mean(LVP_no_offset.es_eivr, 2); mean(LVP_no_offset.eivr_ds, 2);
        mean(LVP_no_offset.ds_ed(1:end-1, :), 2)];
    output.LVP_no_offset_ste = [std(LVP_no_offset.ds_ed(end,:))/sqrt(length(LVP_no_offset.ds_ed(end,:)));
        std(LVP_no_offset.ed_eivc,0,2)/sqrt(length(LVP_no_offset.ed_eivc));
        std(LVP_no_offset.eivc_es,0,2)/sqrt(length(LVP_no_offset.eivc_es));
        std(LVP_no_offset.es_eivr,0,2)/sqrt(length(LVP_no_offset.es_eivr));
        std(LVP_no_offset.eivr_ds,0,2)/sqrt(length(LVP_no_offset.eivr_ds));
        std(LVP_no_offset.ds_ed(1:end-1, :),0,2)/sqrt(length(LVP_no_offset.ds_ed(1:end-1, :)))];
else
    output.LVP_average = [mean(LVP.ds_ed(end,:)); mean(LVP.ed_eivc, 2);
        mean(LVP.eivc_es, 2); mean(LVP.es_ds, 2); mean(LVP.ds_ed(1:end-1, :), 2)];
    output.LVP_ste = [std(LVP.ds_ed(end,:))/sqrt(length(LVP.ds_ed(end,:)));
        std(LVP.ed_eivc,0,2)/sqrt(length(LVP.ed_eivc));
        std(LVP.eivc_es,0,2)/sqrt(length(LVP.eivc_es));
        std(LVP.es_ds,0,2)/sqrt(length(LVP.es_ds));
        std(LVP.ds_ed(1:end-1, :),0,2)/sqrt(length(LVP.ds_ed(1:end-1, :)))];
    output.LVP_no_offset_average = [mean(LVP_no_offset.ds_ed(end,:));
        mean(LVP_no_offset.ed_eivc, 2); mean(LVP_no_offset.eivc_es, 2);
        mean(LVP_no_offset.es_ds, 2); mean(LVP_no_offset.ds_ed(1:end-1, :), 2)];
    output.LVP_no_offset_ste = [std(LVP_no_offset.ds_ed(end,:))/sqrt(length(LVP_no_offset.ds_ed(end,:)));
        std(LVP_no_offset.ed_eivc,0,2)/sqrt(length(LVP_no_offset.ed_eivc));
        std(LVP_no_offset.eivc_es,0,2)/sqrt(length(LVP_no_offset.eivc_es));
        std(LVP_no_offset.es_ds,0,2)/sqrt(length(LVP_no_offset.es_ds));
        std(LVP_no_offset.ds_ed(1:end-1, :),0,2)/sqrt(length(LVP_no_offset.ds_ed(1:end-1, :)))];
end
output.pwp_toggle = haemo.pwp_toggle;

%% Save images of analyses to jpeg files.
disp('Saving analysis figures as JPEG files...');
temp = regexp(mri.volume_path, '[/]', 'split');
output.simulation_name = temp{1};
figure(imageFig);
plot(mri.t, output.LVP_average, 'k*-');
title(['LV trace ', output.simulation_name], 'interpreter', 'none');
xlabel('Time (s)');
ylabel('Pressure (kPa)');
axis([0 1.5 0 20]);
print(imageFig, '-djpeg', ['AnalysisFigures/',output.simulation_name,'_LVP']);
figure(PVFig);
plot(mri.V, output.LVP_average, 'k*-');
title(['LV offsetted PV loop ', output.simulation_name], 'interpreter', 'none');
xlabel('Volume (mL)');
ylabel('Pressure (kPa)');
axis([20 450 0 20]);
print(PVFig, '-djpeg', ['AnalysisFigures/', output.simulation_name, '_PV']);
figure(tracesFig)
if haemo.pwp_toggle == 1
    for i = 1:length(eIVR.t)
        plot([eIVR.t(i), eIVR.t(i)], [0 25], 'k-');
    end
end
title(['Illustration of haemo analyses ', output.simulation_name], 'interpreter', 'none');
xlabel('Time (s)');
ylabel('Pressure (kPa)');
axis([0 1.5 0 25]);
print(tracesFig, '-djpeg', ['AnalysisFigures/', output.simulation_name, '_analyses']);
end