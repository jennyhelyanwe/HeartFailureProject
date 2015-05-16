function AOPressure=ExtractAorticPressure(study_name, debug)

%% This function is designed to extract aortic pressure curves from STF 
%% recordings and overlay them with left ventricular pressures
%%%%%%%%%%%%%%%%%%%%%%%%%% Read in Pressure Data from xls %%%%%%%%%%%%%%%%%
%% Read in the excel data
[Pressure] = xlsread(study_name);

%% Locate all stationary points identified by Jane's (S) and (D)
Stationary_Points=find(isnan(Pressure(:,1)));

if debug
    figure
    plot(Pressure);
    hold on
    plot(Stationary_Points, Pressure(Stationary_Points,2), 'b*');
end

%% Ask user to enter the starting point and number of cycle
SP_Start=input('Please enter the starting stationary point....\n');
No_cycle=input('Please enter the number of cycles....\n');

AP=[];
temporal_spacing_Haemo=4; %% ms
no_sample_points=[];

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(4)/3 scrsz(2)/3 scrsz(3)/3 scrsz(3)])
%% Isolate pressure from multiple cycles
for i=1:No_cycle
    if size(Pressure,2)==3
        AP_Cycle=Pressure(Stationary_Points(SP_Start):(Stationary_Points(SP_Start+2)-1),3);
    else
        AP_Cycle=Pressure(Stationary_Points(SP_Start):(Stationary_Points(SP_Start+2)-1),2);
    end
    no_sample_points_cycle=size(AP_Cycle,1);
    if i>1
        %% Truncate data to keep consistency in overall R-R interval
        if no_sample_points_cycle>no_sample_points(1)
            AP_Cycle((no_sample_points(1)+1):no_sample_points_cycle)=[];
            no_sample_points_cycle=size(AP_Cycle,1);
        %% Append data by calculating the average from the previous cycles
        elseif no_sample_points_cycle<no_sample_points(1)
            AP_Cycle((no_sample_points_cycle+1):no_sample_points(1))=mean(AP(((no_sample_points_cycle+1):no_sample_points(1)),:),2);
            no_sample_points_cycle=size(AP_Cycle,1);
        end
    end
    time=linspace(0,temporal_spacing_Haemo*(no_sample_points_cycle-1),no_sample_points_cycle);
    subplot(3,1,1),plot(time,AP_Cycle,'b*');
    hold on;
    AP=[AP,AP_Cycle];
    SP_Start=SP_Start+2;
    no_sample_points=[no_sample_points;no_sample_points_cycle];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate Average Aortic Pressure %%%%%%%%%%%%%%%%
%% Calculate the average pressure
AP_Average=mean(AP,2);
R_R_Haemo=temporal_spacing_Haemo*(size(AP_Average,1)-1);
TT_Haemo=linspace(0,R_R_Haemo,size(AP_Average,1));
subplot(3,1,1),plot(TT_Haemo,AP_Average,'r+');
hold on;
xlabel('Time (ms)','FontSize',12);
ylabel('Aortic Pressure (mmHg)','FontSize',12);
title('HaemoData Before Temporal Scaling','FontSize',12);
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate dP(ao)/dt %%%%%%%%%%%%%%%%%%%%%%%%%%%
AP_deriv_1=zeros(1,(size(AP_Average,1)));
for i=2:size(AP_deriv_1,2)
    %AP_deriv_1(i)=(2*AP_Average(i)-AP_Average(i+1)-AP_Average(i-1))/8;
    AP_deriv_1(i)=(AP_Average(i)-AP_Average(i-1))/temporal_spacing_Haemo; %forward
    %difference
end
TT_Haemo_deriv_1=linspace(TT_Haemo(1),R_R_Haemo,size(AP_deriv_1,2));
subplot(3,1,2),plot(TT_Haemo_deriv_1,AP_deriv_1,'b*');
grid on;
title(['Aortic Pressure 1st Derivative Study ',study_name(32:37)],'FontSize',12);
xlabel('Time (ms)','FontSize',12);
ylabel('dp/dt','FontSize',12);
hold on;
%{
%% Find the point of deflection which represent (ES and EndIVC)
AP_Average_Sign=AP_deriv_1;
AP_Average_Sign(AP_Average_Sign<0)=0;
AP_Average_Sign(AP_Average_Sign>0)=1;
pt_def=find(diff(AP_Average_Sign)~=0);
%% Plot the point of deflection on the same plot as dP/dt
subplot(3,1,2),plot(TT_Haemo_deriv_1(pt_def(1)),AP_deriv_1(pt_def(1)),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
%subplot(3,1,2),plot(TT_Haemo_deriv_1(pt_def(2)),AP_deriv_1(pt_def(2)),'ko','MarkerSize',8,'MarkerFace','c');
%hold on;
subplot(3,1,2),plot(TT_Haemo_deriv_1(pt_def(3)),AP_deriv_1(pt_def(3)),'ko','MarkerSize',8,'MarkerFace','g');
%% Plot the point of deflection on the same plot as aortic pressure
subplot(3,1,1),plot(TT_Haemo(pt_def(1)),AP_Average(pt_def(1)),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
%subplot(3,1,1),plot(TT_Haemo(pt_def(2)),AP_Average(pt_def(2)),'ko','MarkerSize',8,'MarkerFace','c');
%hold on;
subplot(3,1,1),plot(TT_Haemo(pt_def(3)),AP_Average(pt_def(3)),'ko','MarkerSize',8,'MarkerFace','g');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate dP2(ao)/dt2 %%%%%%%%%%%%%%%%%%%%%%%%%%
AP_deriv_2=zeros(1,(size(AP_Average,1)));
for i=1:(size(AP_deriv_2,2)-1)
   %AP_deriv_2(i)=(2*AP_deriv_1(i)-AP_deriv_1(i+1)-AP_deriv_1(i-1))/8;
   AP_deriv_2(i+1)=(AP_deriv_1(i+1)-AP_deriv_1(i))/temporal_spacing_Haemo; %forward solve
end

%% Apply a smoother to the 2nd derivatives
AP_deriv_2_sm=zeros(1,(size(AP_deriv_1,2)-1));
for i=1:size(AP_deriv_2,2)
    if i>2 & i<(size(AP_deriv_2,2)-2)
        AP_deriv_2_sm(i)=(2*AP_deriv_2(i-2)-AP_deriv_2(i-1)+AP_deriv_2(i+1)+2*AP_deriv_2(i+2))/(0.004*10);
    else
        AP_deriv_2_sm(i)=AP_deriv_2(i);
    end
end
AP_deriv_2=AP_deriv_2_sm;

TT_Haemo_deriv_2=linspace(TT_Haemo(1),R_R_Haemo,size(AP_deriv_2,2));
subplot(3,1,3),plot(TT_Haemo_deriv_2,AP_deriv_2,'k*-');
title(['Aortic Pressure 2nd Derivative Study ', study_name(32:37)],'FontSize',12);
xlabel('Time (ms)','FontSize',12);
ylabel('dp^2/dt^2','FontSize',12);
hold on;
grid on;

%%%%%%%%%%%%%%%%%%%%%% Find AOP at ES and End IVC %%%%%%%%%%%%%%%%%%%%%%%%%
%% Find the minimum AOP 
EndIVC_index=find(AP_Average==min(AP_Average));
%{
%% Find the minimum d2p/dt2 over the first 200ms
AP_deriv_2_Partial=AP_deriv_2(1:50);
ES_index=find(AP_deriv_2_Partial==min(AP_deriv_2_Partial));
%}
%% Find the maximum d2p/dt2 over the first 200ms
AP_deriv_2_Partial=AP_deriv_2(1:50);
ES_index=find(AP_deriv_2_Partial==max(AP_deriv_2_Partial));
pt_def=[EndIVC_index(1),ES_index(1)];
%% Plot the point of deflection on the same plot as aortic pressure
subplot(3,1,1),plot(TT_Haemo(pt_def(1)),AP_Average(pt_def(1)),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
subplot(3,1,2),plot(TT_Haemo_deriv_1(pt_def(1)),AP_deriv_1(pt_def(1)),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
subplot(3,1,3),plot(TT_Haemo_deriv_2(pt_def(1)),AP_deriv_2(pt_def(1)),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
%% Plot the point of deflection on the same plot as aortic pressure
subplot(3,1,1),plot(TT_Haemo(pt_def(2)),AP_Average(pt_def(2)),'ko','MarkerSize',8,'MarkerFace','g');
hold on;
subplot(3,1,2),plot(TT_Haemo_deriv_1(pt_def(2)),AP_deriv_1(pt_def(2)),'ko','MarkerSize',8,'MarkerFace','g');
hold on;
subplot(3,1,3),plot(TT_Haemo_deriv_2(pt_def(2)),AP_deriv_2(pt_def(2)),'ko','MarkerSize',8,'MarkerFace','g');


%%%%%%%%%%%%%%%%%%%%%%%%%%% Output Aortic Pressure %%%%%%%%%%%%%%%%%%%%%%%
AOP_Average_TT(1,:)=TT_Haemo;
AOP_Average_TT(2,:)=AP_Average;
AOPressure.AOP_Average_HaemoTT=AOP_Average_TT;
AOPressure.EndStages(1,:)=[pt_def(1),pt_def(2)];
AOPressure.EndStages(2,:)=[AP_Average(pt_def(1)),AP_Average(pt_def(2))];

return