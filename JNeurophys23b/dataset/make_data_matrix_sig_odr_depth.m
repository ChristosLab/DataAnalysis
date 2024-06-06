% Make dataset of Depth/layer of responsive PFC neurons in ODR task
% J Zhu, 20230608
%% read target neuron numbers from excel
clearvars
excel_name = 'AllNeurons2016FirstCohort';
sheet_name = 'sig';
neuron_wanted = readtable([excel_name '.xlsx'],'sheet',sheet_name);
%% load all neuron dataset
load('odr_data_depth_20230906_raw_max.mat');
file_name = string(neuron_info.Filename);
neuron_info.ID = extract(file_name,1);
%% seg data
[join_T,~,iright] = innerjoin(neuron_wanted,neuron_info);
odr_data = odr_data(iright,:);
%% save dataset
script = char(fread(fopen([mfilename, '.m'])))';
neuron_info = join_T;
save("sig_odr_data_depth_20230906_raw_max.mat","odr_data","neuron_info","script");
disp('Finish running')
clearvars
