% Depth/layer of PFC neurons in ODR task
% find percentage of delay neurons
% J Zhu, 20230831
%% load all neuron data
clearvars
load('odr_data_depth_20230619_raw_max.mat');
file_name = string(neuron_info.Filename);
neuron_info.ID = extract(file_name,1);
%% seg data/label groups
neuron_info.group(neuron_info.Depth<=800)=1;
neuron_info.group(neuron_info.Depth>800&neuron_info.Depth<=1200)=2;
neuron_info.group(neuron_info.Depth>1200)=3;
odr_data = odr_data(~isnan(neuron_info.Depth),:);
neuron_info = neuron_info(~isnan(neuron_info.Depth),:);
%% find numbers
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
nn_all = [];
for g = 1:size(group,1)
    for s = 1:size(stage,1)
        nn_all(g,s) = numel(find(neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)))));
    end
end
%% read target neuron numbers from excel
excel_name = 'neuron_sig_epoch';
sheet_name = 'sig';
neuron_wanted = readtable([excel_name '.xlsx'],'sheet',sheet_name);
% delay_neuron = neuron_wanted(neuron_wanted.Delay==1,:);
%% seg data
[join_T,~,iright] = innerjoin(neuron_wanted,neuron_info,"Keys","Neuron");
%% find numbers
group = unique(join_T.group);
stage = unique(join_T.Stage);
nn_sig = [];
for g = 1:size(group,1)
    for s = 1:size(stage,1)
        nn_sig(g,s) = numel(find(join_T.Neuron(join_T.group==group(g)&contains(join_T.Stage,stage(s)))));
    end
end
%% Percentage of delay neurons
perc = nn_sig./nn_all;