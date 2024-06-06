% neuron firing rate analysis of PFC neurons at diff layers in ODR task
% For ODR task, sig neurons, rate at different epoches
% 20230706, Junda Zhu
%% load data
clearvars
load('sig_odr_data_depth_20230619_raw_max.mat');
file_name = string(neuron_info.Filename);
neuron_info.ID = extract(file_name,1);
%% seg data/label groups
neuron_info.group(neuron_info.Depth<=800)=1;
neuron_info.group(neuron_info.Depth>800&neuron_info.Depth<=1200)=2;
neuron_info.group(neuron_info.Depth>1200)=3;
odr_data = odr_data(~isnan(neuron_info.Depth),:);
neuron_info = neuron_info(~isnan(neuron_info.Depth),:);
%% compute evoked firing rate
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
% epoch_start = 2.0;
% epoch_end= 2.3;
for g = 1:size(group,1)
    nn = [];
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_class = neuron_info.cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        rate_test = [];
        for n = 1:length(Neurons)
            MatData = odr_data_group(n,:);
            rate_sac = [];
            rate_del = [];
            rate_cue = [];
            rate_fix = [];
            for j = 1:length(MatData)
                try
                    rate_sac(j) = mean([MatData{j}.sacrate]);
                    rate_del(j) = mean([MatData{j}.cuedelay]);
                    rate_cue(j) = mean([MatData{j}.cuerate]);
%                     [FR_temp, ntrs_temp] = Get_FRbyneuron_AllTrials_alignCue(MatData,j,epoch_start,epoch_end);
%                     rate_temp(j) = mean(FR_temp);
                    if isfield(MatData{j},'fixrate')
                        rate_fix(j) = mean([MatData{j}.fixrate]);
                    else
                        rate_fix(j) = mean([MatData{j}.fix]);
                    end
                catch
                    lasterr
                    rate_sac(j) = nan;
                    rate_del(j) = nan;
                    rate_cue(j) = nan;
                    rate_fix(j) = nan;
                end
            end
%                         rate_test(n) = mean(rate_fix,"omitnan"); % from all trial (for baseline)
            rate_test(n) = rate_del(Best_class(n))-rate_fix(Best_class(n));
        end
        rate_test_mean = mean(rate_test,'omitnan');
        rate_test_sem = std(rate_test,'omitnan')./sqrt(length(rate_test));
        nn = size(odr_data_group,1);
        % save for plot
        plt_save.rate{g,s} = rate_test;
        plt_save.rate_sem{g,s} = rate_test_sem;
        plt_save.nn{g,s} = nn;
        plt_save.rate_mean{g,s} = rate_test_mean;
    end
end
%% permutation test
p_value = [];
for g = 1:size(group,1)
    % Prepare data
    data1 = plt_save.rate{g,1};
    data2 = plt_save.rate{g,2};
    % Compute the observed difference
    obs_diff = mean(data2) - mean(data1);
    % Concatenate the data
    all_data = [data1, data2];
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
    p_value(g) = sum(abs(perm_diffs) >= abs(obs_diff)) / num_permutations;
end