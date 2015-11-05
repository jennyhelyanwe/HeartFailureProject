function Align_LVP_AOP

clear all;
clc;
close all;
debug = true;
%
fprintf('\n************ Start to analyse aortic pressure *************\n');
temporal_spacing_Haemo=4;   % 4 ms between each haemodynamic measurement. 

%% Enter the excel filename to be analysed (Aortic Pressure)
study_name='Aortic_Pressure_Separate_Study/STF-19-MR-276276.xlsx';
AOPressure=ExtractAorticPressure(study_name, debug);
% print the interval between peaks for aortic pressure trace
fprintf('===== The duration between peak aortic pressure is %f ======\n',AOPressure.AOP_Average_HaemoTT(1,end));
fprintf('\n***********************************************************\n');

fprintf('\n************ Start to analyse LV pressure *************\n');
%}

%% Enter the excel filename to be analysed (LV Pressure)
study_name='LV_Pressure_Separate_Study/STF-19-MR-276276.xlsx';
LVPressure=ExtractLVPressure(study_name, debug);
fprintf('===== The duration between peak LV pressure is %f ======\n',LVPressure.LVP_Average_HaemoTT(1,end));
fprintf('\n***********************************************************\n');

%% Plot the two traces together
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(4)/3 scrsz(2)/3 scrsz(3)/2.5 scrsz(3)])
subplot(3,1,1),plot(AOPressure.AOP_Average_HaemoTT(1,:),AOPressure.AOP_Average_HaemoTT(2,:),'r*');
hold on;
subplot(3,1,1),plot(LVPressure.LVP_Average_HaemoTT(1,:),LVPressure.LVP_Average_HaemoTT(2,:),'b*');
grid on;
subplot(3,1,1),plot(AOPressure.AOP_Average_HaemoTT(1,AOPressure.EndStages(1,1)),AOPressure.AOP_Average_HaemoTT(2,AOPressure.EndStages(1,1)),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
subplot(3,1,1),plot(AOPressure.AOP_Average_HaemoTT(1,AOPressure.EndStages(1,2)),AOPressure.AOP_Average_HaemoTT(2,AOPressure.EndStages(1,2)),'ko','MarkerSize',8,'MarkerFace','g');
hold on;

xlabel('Time (ms)','FontSize',12);
ylabel('Pressure (mmHg)','FontSize',12);
title(['LV & Aortic Pressure Study ',study_name(28:33)],'FontSize',12);
legend('Aortic Pressure','LV Pressure','EndIVC','ES');

%% Select LVP at EndIVC, ES and EndIVR
% EndIVC
Pressure_diff=abs(LVPressure.LVP_Average_HaemoTT(2,:)-AOPressure.EndStages(2,1));
LVP_EndIVC_index=find(Pressure_diff==min(Pressure_diff(size(Pressure_diff,2)/2:size(Pressure_diff,2))));
LVP_EndIVC=LVPressure.LVP_Average_HaemoTT(2,LVP_EndIVC_index);
% ES
Pressure_diff=abs(LVPressure.LVP_Average_HaemoTT(2,:)-AOPressure.EndStages(2,2));
% This limit the search to the 1st part of the data only
LVP_ES_index=find(Pressure_diff==min(Pressure_diff(1:floor(size(Pressure_diff,2)/2))));
LVP_ES=LVPressure.LVP_Average_HaemoTT(2,LVP_ES_index);
% Collect LVP at all end stages
LVPressure.EndStages(1,3)=LVP_EndIVC_index;
LVPressure.EndStages(1,4)=LVP_ES_index;
LVPressure.EndStages(2,3)=LVP_EndIVC;
LVPressure.EndStages(2,4)=LVP_ES;
LVPressure.EndStages(3,:)=(LVPressure.EndStages(1,:)-1).*temporal_spacing_Haemo;

subplot(3,1,1),plot(LVPressure.LVP_Average_HaemoTT(1,LVPressure.EndStages(1,1)),LVPressure.LVP_Average_HaemoTT(2,LVPressure.EndStages(1,1)),'ko','MarkerSize',8,'MarkerFace','m');
hold on;
subplot(3,1,1),plot(LVPressure.LVP_Average_HaemoTT(1,LVPressure.EndStages(1,2)),LVPressure.LVP_Average_HaemoTT(2,LVPressure.EndStages(1,2)),'ko','MarkerSize',8,'MarkerFace','c');
hold on;
subplot(3,1,1),plot(LVPressure.LVP_Average_HaemoTT(1,LVPressure.EndStages(1,3)),LVPressure.LVP_Average_HaemoTT(2,LVPressure.EndStages(1,3)),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
subplot(3,1,1),plot(LVPressure.LVP_Average_HaemoTT(1,LVPressure.EndStages(1,4)),LVPressure.LVP_Average_HaemoTT(2,LVPressure.EndStages(1,4)),'ko','MarkerSize',8,'MarkerFace','g');
hold on;
subplot(3,1,1),plot(LVPressure.LVP_Average_HaemoTT(1,LVPressure.EndStages(1,5)),LVPressure.LVP_Average_HaemoTT(2,LVPressure.EndStages(1,5)),'ko','MarkerSize',8,'MarkerFace','k');
hold on;
%legend('Aortic Pressure','LV Pressure','EndIVC-AOP','ES-AOP','DS-LVP','ED-LVP','EndIVC-LVP','ES-LVP','EndIVR-LVP','Location','EastOutside','Orientation','vertical');
grid on;

%% Shift the end stages to start from ED
LVPressure.EndStages=circshift(LVPressure.EndStages',-1)';

%%%%%%%%% Temporal alignment of the traces %%%%%%%%%%%%%%%%%%%%%%
%{
%% Temporal scale the aortric trace
AOP_LVP_ratio=LVPressure.LVP_Average_HaemoTT(1,end)\AOPressure.AOP_Average_HaemoTT(1,end);
AOPressure.AOP_Average_HaemoTT(3,:)=AOPressure.AOP_Average_HaemoTT(1,:).*AOP_LVP_ratio;
for i=1:size(LVPressure.LVP_Average_HaemoTT(1,:),2)
    TT_diff=abs(AOPressure.AOP_Average_HaemoTT(3,:)-LVPressure.LVP_Average_HaemoTT(1,i));
    match_index=find(TT_diff==min(TT_diff));
    AOPressure.AOP_Average_HaemoTT(4,i)=AOPressure.AOP_Average_HaemoTT(2,match_index);
end
%% Remove the zero components
non_zero_index=find(AOPressure.AOP_Average_HaemoTT(4,:)~=0);
AOPressure.AOP_Average_Scaled=AOPressure.AOP_Average_HaemoTT(4,non_zero_index);
plot(LVPressure.LVP_Average_HaemoTT(1,:),AOPressure.AOP_Average_Scaled,'c*');
%}
MRI_file='LV_Pressure_Separate_Study/STF-19_MRI_info.mat';
study_name='STF_19';
ExtractCavityVolume(MRI_file,LVPressure,study_name);

return