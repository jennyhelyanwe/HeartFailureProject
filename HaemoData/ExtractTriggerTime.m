function [RR_mean]=ExtractTriggerTime(image_name,no_frames,current_study)


%% This function is designed to extract trigger time from all dicom images

%% Change into the image directory
cd(image_name);
all_files=dir;
all_files(1:2,:)=[];

%% Calculate number of slices
no_slices=size(all_files,1)/no_frames;
k=1;
TT_all_frames=zeros(no_frames,no_slices);
TS=zeros(no_slices,1);

for i=1:no_slices
    for j=1:no_frames
        infor=dicominfo(all_files(k).name);
        TT_all_frames(j,i)=infor.TriggerTime;
        k=k+1;
    end
    %% Sort the trigger time to see whether there is image out of place
    [TT_sort,index]=sortrows(TT_all_frames(:,i));
    frame_number=linspace(1,no_frames,no_frames);
    frame_diff=index-frame_number';
    if ~isempty(find(frame_diff~=0))
        fprintf('*** Warning: Slice %d requires reordering of the frames ......\n',i);
        TT_all_frames(:,i)=TT_sort;
    end
    
    %% Calculate the temporal spacing
    TS_all_frames=diff(TT_sort);
    if diff(TS_all_frames)<1e-4
        TS(i)=TS_all_frames(1);
        fprintf('*** Slice %d has consistent temporal spacing ......\n',i);
    else
        fprintf('*** Warning: Slice %d does not have consistent temporal spacing ......\n',i);
    end
end

%% Calculate the mean temporal spacing
TS_mean=mean(TS);
TS_std=std(TS);
RR_mean=TS_mean*(no_frames-1);

fprintf('+++++ The mean temporal spacing is %f .....\n',TS_mean);
fprintf('+++++ The standard deviation for temporal spacing is %f .....\n',TS_std);
fprintf('+++++ The mean R-R interval is %f .........\n',RR_mean);

MRI_Info.RR_mean=RR_mean;
MRI_Info.TS_mean=TS_mean;
study_file=['../../HaemoDataFromSTF/Pressure_Separate_Study/',current_study,'_MRI_info'];
save(study_file,'MRI_Info'); 

cd('../../HaemoDataFromSTF');

return
