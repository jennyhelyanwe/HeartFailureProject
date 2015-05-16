function ExtractCavityVolume(ImageFile,LVP,study_name)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Read in R-R Interval for MRI %%%%%%%%%%%%%%%%%%
%% Read in the R-R interval for the MRI data
isexist=exist(ImageFile);
if isexist==0
    cd('STF_Images/');
    image_name=input('Please enter the name of the image study .....\n','s');
    No_Frame=input('Please enter the total number of frames for the current study .....\n');
    current_study=ImageFile(28:33);
    [R_R_MRI]=ExtractTriggerTime(image_name,No_Frame,current_study);
    TT_MRI=R_R_MRI/(No_Frame-1);
    fprintf('+++++ The R-R interval is %f and the timing between each frame is %f.\n',R_R_MRI, TT_MRI);
else
    MRI_Info=load(ImageFile);
    R_R_MRI=MRI_Info.MRI_Info.RR_mean;
    TT_MRI=MRI_Info.MRI_Info.TS_mean;
    fprintf('+++++ The mean R-R interval is %f, and TT interval is %f.........\n',R_R_MRI, TT_MRI);
end

%{
fprintf('***** The mean R-R interval for pressure is %f .....\n',R_R_Haemo);
%%%%%%%%%% Temporally scale the Haemo data to match with MRI %%%%%%%%%%%%%%
%% Scale the R-R interval for the Haemo data
R_R_Haemo_Scale=R_R_MRI/R_R_Haemo;
%% Sacle the triger time for the Haemo data
LVP_Average_TT(3,:)=LVP_Average_TT(1,:).*R_R_Haemo_Scale;
%}
%}

%% At this stage, LVP.LVP_Average_HaemoTT is from peak-peak
%% Reorder the pressure tracet to start from ED to ED
no_sample=size(LVP.LVP_Average_HaemoTT,2);
LVP_ED_EndCycle=LVP.LVP_Average_HaemoTT(:,LVP.EndStages(1,1):no_sample);
LVP_Start_ED=LVP.LVP_Average_HaemoTT(:,1:(LVP.EndStages(1,1)-1));
LVP.LVP_Average_HaemoTT_Reorder=[LVP_ED_EndCycle,LVP_Start_ED]; %% ED-ED
%% Find the index for the end stage pressure
for k=1:size(LVP.EndStages,2)
    Index_reorder=find(LVP.LVP_Average_HaemoTT_Reorder(2,:)==LVP.EndStages(2,k));
    LVP.EndStages_Reorder(k)=Index_reorder(1);
end

%%%%%%%%%%% Locate pressure at the end of each cardiac phases %%%%%%%%%%%%
No_Frame=input('Please enter the total number of frames for the current study .....\n');
%% Locate the time for ED, End_IVC, ES, End_IVR and DS
ED_Frame=input('Please enter the frame for the ED phase .....\n');
ED=(ED_Frame-1)*TT_MRI;
End_IVC_Frame=input('Please enter the frame for the End IVC phase .....\n');
End_IVC=(End_IVC_Frame-1)*TT_MRI;
ES_Frame=input('Please enter the frame for the ES phase .....\n');
ES=(ES_Frame-1)*TT_MRI;
End_IVR_Frame=input('Please enter the frame for the End IVR phase .....\n');
End_IVR=(End_IVR_Frame-1)*TT_MRI;
DS_Frame=input('Please enter the frame for the DS phase .....\n');
DS=(DS_Frame-1)*TT_MRI;

LVP.EndStages(4,:)=[1,End_IVC_Frame,ES_Frame,End_IVR_Frame,DS_Frame];

no_sample=size(LVP.LVP_Average_HaemoTT,2);
LVP_MRI_Matched=[];
Sample_Previous=1;

for p=1:size(LVP.EndStages,2)
    if p<5
        %% Interpolate the pressure between stages
        LVP_Inter_Frame=LVP.EndStages(4,p+1)-LVP.EndStages(4,p);
        LVP_Inter_Sample=LVP.EndStages_Reorder(1,p+1)-LVP.EndStages_Reorder(1,p);
        LVP_Inter=LVP.LVP_Average_HaemoTT_Reorder(:,LVP.EndStages_Reorder(1,p):LVP.EndStages_Reorder(1,p+1));
    else
        LVP_Inter_Frame=No_Frame-LVP.EndStages(4,p);
        LVP_Inter_Sample=no_sample-LVP.EndStages_Reorder(1,p);
        LVP_Inter=LVP.LVP_Average_HaemoTT_Reorder(:,LVP.EndStages_Reorder(1,p):no_sample);
    end
   
    LVP_temp_step=LVP_Inter_Sample/LVP_Inter_Frame;
    LVP_Inter_Phase=zeros(1,LVP_Inter_Frame+1);
    LVP_Inter_Phase(1,1)=Sample_Previous;
    LVP_Inter_Phase(1,end)=Sample_Previous+LVP_Inter_Sample;
    LVP_Inter_Phase(2,1)=LVP_Inter(2,1);
    LVP_Inter_Phase(2,end)=LVP_Inter(2,end);

    for k=2:LVP_Inter_Frame
        LVP_Inter_Phase(1,k)=Sample_Previous+floor((k-1)*LVP_temp_step);
        LVP_Inter_Phase(2,k)=LVP_Inter(2,floor((k-1)*LVP_temp_step)); 
    end
    
    Sample_Previous=LVP_Inter_Phase(1,end);  

    if p==1
        LVP_MRI_Matched=[LVP_MRI_Matched,LVP_Inter_Phase];
    else
        LVP_MRI_Matched=[LVP_MRI_Matched,LVP_Inter_Phase(:,2:end)];
    end
end


%% Plot the interpolated LVP
subplot(3,1,2),plot(LVP.LVP_Average_HaemoTT(1,:),LVP.LVP_Average_HaemoTT_Reorder(2,:),'b*','MarkerSize',8);
hold on;
LVP.LVP_Average_Resampled(1,:)=LVP.LVP_Average_HaemoTT(1,LVP_MRI_Matched(1,:));
LVP.LVP_Average_Resampled(2,:)=LVP_MRI_Matched(2,:);

%% Check for peak discontinuity in interpolated LVP
peak_index = find(LVP.LVP_Average_Resampled(2,:) == max(LVP.LVP_Average_Resampled(2,:)));
if (abs(LVP.LVP_Average_Resampled(2,peak_index) -LVP.LVP_Average_Resampled(2,peak_index-1)) > 7)
    pre_discont = peak_index -1;
    post_discont = peak_index;
    %LVP.LVP_Average_Resampled(2,:) = Peak_Interpolation(LVP.LVP_Average_Resampled(2,:), pre_discont, post_discont);
elseif (abs(LVP.LVP_Average_Resampled(2,peak_index) -LVP.LVP_Average_Resampled(2,peak_index+1)) > 7)
    pre_discont = peak_index;
    post_discont = peak_index + 1;
    %LVP.LVP_Average_Resampled(2,:) = Peak_Interpolation(LVP.LVP_Average_Resampled(2,:), pre_discont, post_discont);
end

subplot(3,1,2),plot(LVP.LVP_Average_Resampled(1,:),LVP.LVP_Average_Resampled(2,:),'ro','MarkerSize',6,'MarkerFace','w');
hold on;
grid on;
%legend('LVP-Average','LVP@MRI Frame','Location','SouthOutside','Orientation','horizontal');
legend('LVP-Average','LVP@MRI Frame');


%% Convert pressure to kPa
LVP.LVP_Average_Resampled(2,:)=LVP_MRI_Matched(2,:)*101.32/760;
%% Find the minimum pressure (DS pressure) and offset the pressure
LVP.LVP_Average_Resampled(3,:)=LVP.LVP_Average_Resampled(2,:)-min(LVP.LVP_Average_Resampled(2,:));

%% Deal with discontinuity
%LVP.LVP_Average_Resampled(3,:) = Peak_Interpolation(LVP.LVP_Average_Resampled(3,:), pre_discont, post_discont);
%% Write the register pressure to a text file
filenamew=(['RegisteredPressure/',study_name,'_registered_LVP.txt']);
fidw=fopen(filenamew,'w');
for i=1:size(LVP.LVP_Average_Resampled(2,:),2)
    fprintf(fidw,'%d    %f  \n',i,LVP.LVP_Average_Resampled(3,i));
end
fclose(fidw);

%%%%%%%%%%%%%% Plot pressure volume curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CIM_Volume_File='C:\CIMSoftware\CIMDATA\STF_016\volumes\info\STF_016_model_ca.model_STF_016_AY';
CIM_Volume_File=(['../CIM_Models/Studies/',study_name,'/volumes/info/',study_name,'_model_ca.model_',study_name,'_AY'])
fidr=fopen(CIM_Volume_File,'r');
%% Initialise the line being read
line_read=0;
junk=fgets(fidr);
line_read=line_read+1;
while isempty(strfind(junk, 'TOTAL'))
    junk=fgets(fidr);
    line_read=line_read+1;
end
junk = fgets(fidr);
junk = fgets(fidr);
volume_name=fscanf(fidr,'%s',2);
CIM_GPM_Volume=fscanf(fidr,'%f',No_Frame);
fclose(fidr);
%%% Shift the cavity volume so that volume starts from ED
%if ED_Frame>End_IVC_Frame
%    CIM_GPM_Volume=circshift(CIM_GPM_Volume,(No_Frame-ED_Frame)+1);
%end

%{
%% Plot pressure-volume loop
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(4)/3 scrsz(2)/3 scrsz(3)/3 scrsz(3)])

subplot(3,1,1),plot(LVP_Average_Frame(1,:),LVP_Average_Frame(2,:),'b*');
hold on;
subplot(3,1,1),plot(LVP_Average_Frame(1,ED_Frame),LVP_Average_Frame(2,ED_Frame),'ko','MarkerSize',8,'MarkerFace','r');
hold on;
subplot(3,1,1),plot(LVP_Average_Frame(1,End_IVC_Frame),LVP_Average_Frame(2,End_IVC_Frame),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
subplot(3,1,1),plot(LVP_Average_Frame(1,ES_Frame),LVP_Average_Frame(2,ES_Frame),'ko','MarkerSize',8,'MarkerFace','g');
hold on;
subplot(3,1,1),plot(LVP_Average_Frame(1,End_IVR_Frame),LVP_Average_Frame(2,End_IVR_Frame),'ko','MarkerSize',8,'MarkerFace','b');
hold on;
subplot(3,1,1),plot(LVP_Average_Frame(1,DS_Frame),LVP_Average_Frame(2,DS_Frame),'ko','MarkerSize',8,'MarkerFace','c');
title('LV Pressure Temporal Profile Study STF-9 ','FontSize',12);
xlabel('Frame','FontSize',12);
ylabel('LV Pressure (kPa)','FontSize',12);
grid on;

subplot(3,1,2),plot(LVP_Average_Frame(1,:),CIM_GPM_Volume,'b*');
hold on;
subplot(3,1,2),plot(LVP_Average_Frame(1,ED_Frame),CIM_GPM_Volume(ED_Frame),'ko','MarkerSize',8,'MarkerFace','r');
hold on;
subplot(3,1,2),plot(LVP_Average_Frame(1,End_IVC_Frame),CIM_GPM_Volume(End_IVC_Frame),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
subplot(3,1,2),plot(LVP_Average_Frame(1,ES_Frame),CIM_GPM_Volume(ES_Frame),'ko','MarkerSize',8,'MarkerFace','g');
hold on;
subplot(3,1,2),plot(LVP_Average_Frame(1,End_IVR_Frame),CIM_GPM_Volume(End_IVR_Frame),'ko','MarkerSize',8,'MarkerFace','b');
hold on;
subplot(3,1,2),plot(LVP_Average_Frame(1,DS_Frame),CIM_GPM_Volume(DS_Frame),'ko','MarkerSize',8,'MarkerFace','c');
title('LV Cavity Volume Temporal Profile Study STF-9 ','FontSize',12);
xlabel('Frame','FontSize',12);
ylabel('LV Cavity Volume (ml)','FontSize',12);
grid on;

%}
subplot(3,1,3),plot(CIM_GPM_Volume',LVP.LVP_Average_Resampled(2,:)','b*');
%axis([min(CIM_GPM_Volume)*0.5 max(CIM_GPM_Volume)*1.2 0 max(LVP.LVP_Average_Resampled(2,:))*1.2]);
axis([0 400 0 max(LVP.LVP_Average_Resampled(2,:))*1.2]);
t = title(['Pressure Volume Loop ',study_name],'FontSize',12);
set(t, 'Interpreter', 'none');
xlabel('LV Cavity Volume (ml)','FontSize',12);
ylabel('LV Pressure (kPa)','FontSize',12);
hold on;
Pressure_Volume=[CIM_GPM_Volume,LVP.LVP_Average_Resampled(2,:)'];
subplot(3,1,3),plot(Pressure_Volume(ED_Frame,1),Pressure_Volume(ED_Frame,2),'ko','MarkerSize',8,'MarkerFace','c');
hold on;
subplot(3,1,3),plot(Pressure_Volume(End_IVC_Frame,1),Pressure_Volume(End_IVC_Frame,2),'ko','MarkerSize',8,'MarkerFace','y');
hold on;
subplot(3,1,3),plot(Pressure_Volume(ES_Frame,1),Pressure_Volume(ES_Frame,2),'ko','MarkerSize',8,'MarkerFace','g');
hold on;
subplot(3,1,3),plot(Pressure_Volume(End_IVR_Frame,1),Pressure_Volume(End_IVR_Frame,2),'ko','MarkerSize',8,'MarkerFace','k');
hold on;
subplot(3,1,3),plot(Pressure_Volume(DS_Frame,1),Pressure_Volume(DS_Frame,2),'ko','MarkerSize',8,'MarkerFace','m');
grid on;

LVP_ED=Pressure_Volume(ED_Frame,2)-Pressure_Volume(DS_Frame,2);
LVP_End_IVC=Pressure_Volume(End_IVC_Frame,2)-Pressure_Volume(DS_Frame,2);
LVP_ES=Pressure_Volume(ES_Frame,2)-Pressure_Volume(DS_Frame,2);
LVP_End_IVR=Pressure_Volume(End_IVR_Frame,2)-Pressure_Volume(DS_Frame,2);
LVP_DS=Pressure_Volume(DS_Frame,2)-Pressure_Volume(DS_Frame,2);

LVP.Pressure_EndStage_AfterOffset=[LVP_ED,LVP_End_IVC,LVP_ES,LVP_End_IVR,LVP_DS];
LVP.PressureVolume=Pressure_Volume;

%save Pressure_Separate_Study/LVP_STF-9 LVP
%}


return