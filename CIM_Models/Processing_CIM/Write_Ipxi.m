function Write_Ipxi

clc;
clear all;

%% This file is designed to write ipxi file

%% Define the file name to be written
filenamew1=('Surface_Points_Endo.ipxi');
fidw1=fopen(filenamew1,'w');

%% Define the total number of data points to be generated
no_elem=16;
no_points_xi1=5;
no_points_xi2=5;
no_data=1;

%% Endocardial surface points
for i=1:4 %% loop through elements 1-4, few number of points for the apical elements
    for j=1:(no_points_xi2)  %% loop through direction xi_2
        for k=1:(no_points_xi1)  %% loop through direction xi_1
            fprintf(fidw1,'  %d     %d     %1.12f     %1.12f     %1.12f     \n',no_data,i,(1/no_points_xi1)*(k-1),(1/no_points_xi2)*(j-1),0.00000);
            no_data=no_data+1;
        end
    end
end

%% Increase points for the non-apical elements
no_points_xi1=9;
no_points_xi2=9;

for i=5:no_elem %% loop through all elements
    for j=1:(no_points_xi2+1)  %% loop through direction xi_2
        for k=1:(no_points_xi1+1)  %% loop through direction xi_1
            fprintf(fidw1,'  %d     %d     %1.12f     %1.12f     %1.12f     \n',no_data,i,(1/no_points_xi1)*(k-1),(1/no_points_xi2)*(j-1),0.00000);
            no_data=no_data+1;
        end
    end
end
fclose(fidw1);

%% Define the file name to be written
filenamew1=('Surface_Points_Epi.ipxi');
fidw1=fopen(filenamew1,'w');

no_points_xi1=5;
no_points_xi2=5;
no_data=1;

%% Epicardial surface points
%% Endocardial surface points
for i=1:4 %% loop through elements 1-4, few number of points for the apical elements
    for j=1:(no_points_xi2)  %% loop through direction xi_2
        for k=1:(no_points_xi1)  %% loop through direction xi_1
            fprintf(fidw1,'  %d     %d     %1.12f     %1.12f     %1.12f     \n',no_data,i,(1/no_points_xi1)*(k-1),(1/no_points_xi2)*(j-1),1.00000);
            no_data=no_data+1;
        end
    end
end

%% Increase points for the non-apical elements
no_points_xi1=9;
no_points_xi2=9;


for i=5:no_elem  %% loop through all elements
    for j=1:(no_points_xi2+1)  %% loop through direction xi_2
        for k=1:(no_points_xi1+1)  %% loop through direction xi_1
            fprintf(fidw1,'  %d     %d     %1.12f     %1.12f     %1.12f   \n',no_data,i,(1/no_points_xi1)*(k-1),(1/no_points_xi2)*(j-1),1.00000);
            no_data=no_data+1;
        end
    end
end


fprintf('Finished writing the ipxi file ......\n');

return