% Make dataset of sulcus unit or not of responsive PFC neurons in ODR task
% J Zhu, 2023828
%% read target neuron numbers from excel
clearvars
excel_name = 'neuron_sulcus_cood';
sheet_name = 'neuron_sulcus_cood';
neuron_wanted = readtable([excel_name '.xlsx'],'sheet',sheet_name);
%% load sig neuron dataset
load('sig_odr_data_depth_20230906_raw_max.mat');
file_name = string(neuron_info.Filename);
neuron_info.ID = extract(file_name,1);
%% add sulcus column
[join_T,~,iright] = innerjoin(neuron_wanted,neuron_info);
neuron_info.sulcus(:) = 0;
neuron_info.sulcus(iright) = 1;
%% save dataset
script = char(fread(fopen([mfilename, '.m'])))';
save("sig_odr_data_depth_20230906_raw_max_sulcus.mat","odr_data","neuron_info","script");
disp('Finish running')
clearvars