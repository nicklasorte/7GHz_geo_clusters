clear;
clc;
close all force;
close all;
app=NaN(1);  %%%%%%%%%This is to allow for Matlab Application integration.
format shortG
%format longG
top_start_clock=clock;
folder1='C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\7GHz Clusters';
cd(folder1)
addpath(folder1)
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\Basic_Functions')
addpath('C:\Local Matlab Data\Local MAT Data') %%%%%%%One Drive Error with mat files
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%Clusters for Kim
tf_repull_p2p=0%1%0%1%0
excel_filename_p2p='iQlinkDB_Report_Link_LIP 2026-2-8_Xlinks.xlsx' 
mat_filename_str_p2p=strcat('7ghz_xlinks_p2p_data.mat')

tic;
[cell_p2p_data]=load_full_excel_rev1(app,mat_filename_str_p2p,excel_filename_p2p,tf_repull_p2p);
toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Pull the name and
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end points
cell_header=cell_p2p_data(1,:)';
col_ulink_idx=find(matches(cell_header,'ULINK_ID'));
col_alat_idx=find(matches(cell_header,'A_LAT'));
col_alon_idx=find(matches(cell_header,'A_LON'));
col_blat_idx=find(matches(cell_header,'B_LAT'));
col_blon_idx=find(matches(cell_header,'B_LON'));


cell_data=cell_p2p_data(2:end,[col_ulink_idx,col_alat_idx,col_alon_idx,col_blat_idx,col_blon_idx]);

% x_km=100
% n_cluster=10 
% cluster_cap=300


%%%%%%%%%%This seems to be the best inputs
x_km=100
n_cluster=8 
cluster_cap=400
iterations=10


%%%%%%%%%%%%%%%%%%This leaves a little florida cluster
% % % x_km=150
% % % n_cluster=8 
% % % cluster_cap=350
% % % iterations=10


% %%%%%%%%%%
% x_km=150
% n_cluster=8 
% cluster_cap=375
% iterations=10



%%%%%%%%Getting close.
outFinal = buildFinalClusters_COG(cell_data, x_km, n_cluster, cluster_cap, iterations);

outFinal.info
outFinal.clusterSizes
%%%%%%%%%%%%%%%%%%%%%%%%%%Plot
plotFinalClustersCONUS(outFinal,strcat('clusters_',num2str(x_km),'km_',num2str(n_cluster),'_',num2str(cluster_cap),'_',num2str(iterations),'.jpg'), 300)

%%%%%%%%%%Save to Excel
saveClustersToExcel(cell_p2p_data, outFinal, strcat('clusters_',num2str(x_km),'km_',num2str(n_cluster),'_',num2str(cluster_cap),'_',num2str(iterations),'.xlsx'));






end_clock=clock;
total_clock=end_clock-top_start_clock;
total_seconds=total_clock(6)+total_clock(5)*60+total_clock(4)*3600+total_clock(3)*86400;
total_mins=total_seconds/60;
total_hours=total_mins/60;
if total_hours>1
    strcat('Total Hours:',num2str(total_hours))
elseif total_mins>1
    strcat('Total Minutes:',num2str(total_mins))
else
    strcat('Total Seconds:',num2str(total_seconds))
end
close all force;
cd(folder1)
'Done'













