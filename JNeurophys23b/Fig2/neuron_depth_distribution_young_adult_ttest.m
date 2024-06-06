% Depth of PFC neurons recorded from young and adult stage. compare
% ditribution of the two stages
% J Zhu, 20230629
%% load data
clearvars
load('sig_odr_data_depth_20230906_raw_max_sulcus.mat');
file_name = string(neuron_info.Filename);
neuron_info.ID = extract(file_name,1);
%% select neuron: optional
odr_data = odr_data(neuron_info.sulcus==0,:);
neuron_info = neuron_info(neuron_info.sulcus==0,:);
%% seg data/label groups
neuron_info.group(neuron_info.Depth<=800)=1;
neuron_info.group(neuron_info.Depth>800&neuron_info.Depth<=1200)=2;
neuron_info.group(neuron_info.Depth>1200)=3;
odr_data = odr_data(~isnan(neuron_info.Depth),:);
neuron_info = neuron_info(~isnan(neuron_info.Depth),:);
id = ['i','j','k','l'];
stage = unique(neuron_info.Stage);
% depth_s = table;
%% Distribution of neurons depth
depth_plt1 = neuron_info.Depth;
group_plt1 = neuron_info.ID;
stage_plt1 = neuron_info.Stage;
for s = 1:size(stage,1)
    depth_s{s} = depth_plt1(contains(stage_plt1,stage(s)));
end
data1 = depth_s{1};
data2 = depth_s{2};
%% ttest
[h,p,ci,stats] = ttest2(depth_s{1},depth_s{2});
%% permutation
% Compute the observed difference
obs_diff = mean(data2) - mean(data1);
% Concatenate the data
all_data = [data1; data2];
% Perform the permutation test
n = length(all_data);
n1 = length(data1);
n2 = length(data2);
num_permutations = 1000;
perm_diffs = zeros(1, num_permutations);
for i = 1:num_permutations
    perm_data = all_data(randperm(n));
    perm_diffs(i) = mean(perm_data(1:n1)) - mean(perm_data((n1+1):end));
end
% Compute the p-value
p_value = sum(abs(perm_diffs) >= abs(obs_diff)) / num_permutations;