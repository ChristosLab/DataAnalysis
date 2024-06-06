% Time-resolved ROC analysis of PFC neurons at diff layers in ODR task
% J Zhu, 20230707
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
%% oppo locations
opp_index = [5 6 7 8 1 2 3 4];
neuron_info.opp_del= opp_index(neuron_info.del)';
neuron_info.opp_cue = opp_index(neuron_info.cue)';
neuron_info.opp_sac= opp_index(neuron_info.sac)';
%% compute firing rate ROC on best cue location/opp location
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
for g = 1:size(group,1)
    nn = [];
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_class = neuron_info.cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Opp_class = neuron_info.opp_cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        neuron_ROC = [];
        rate_best = [];
        rate_opp = [];
        for n = 1:length(Neurons)
            MatData = odr_data_group(n,:);
            rate_sac = [];
            rate_del = [];
            rate_cue = [];
            rate_fix = [];
            for j = 1:length(MatData)
                try
                    rate_sac{j} = [MatData{j}.sacrate];
                    rate_del{j} = [MatData{j}.cuedelay];
                    rate_cue{j} = [MatData{j}.cuerate];
                catch
                    lasterr
                    rate_sac(j) = nan;
                    rate_del(j) = nan;
                    rate_cue(j) = nan;
                end
            end
           neuron_ROC(n) = arrayROC(rate_del{Best_class(n)},rate_del{Opp_class(n)});
        end
        roc_mean = mean(neuron_ROC,'omitnan');
        roc_sem = std(neuron_ROC,'omitnan')./sqrt(length(neuron_ROC));
        nn = size(odr_data_group,1);
        % save for plot
        plt_save.roc{g,s} = neuron_ROC;
        plt_save.roc_sem{g,s} = roc_sem;
        plt_save.nn{g,s} = nn;
        plt_save.roc_mean{g,s} = roc_mean;
    end
end
%% permutation test on firing rate auROC
p_value = [];
for g = 1:size(group,1)
    % Prepare data
    data1 = plt_save.roc{g,1};
    data2 = plt_save.roc{g,2};
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
%% compute time resolved ROC on best cue location/opp location
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
wid = 200; % ms
noedge = 1;
for g = 1:size(group,1)
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_class = neuron_info.cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Opp_class = neuron_info.opp_cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        neuron_ROC = [];
        for n = 1:length(Neurons)
            try
                rate_cue = [];
                rate_opp_cue = [];
                [spiketrain_temp1, ntrs_temp1] = Get_spiketrain_partial_aligncue(odr_data_group(n,:),Best_class(n),[]);
                [spiketrain_temp2, ntrs_temp2] = Get_spiketrain_partial_aligncue(odr_data_group(n,:),Opp_class(n),[]);
                for ntr = 1:min(ntrs_temp1,ntrs_temp2)
                    [spks1, tlo1, thi1] = spkmtx(spiketrain_temp1(ntr),0,[-1000,3500]);
                    t1 = tlo1:thi1;
                    rate1 = 1000*smooth_es(spks1', wid, noedge)';
                    rate_cue(ntr,:) = rate1;
                    [spks2, tlo2, thi2] = spkmtx(spiketrain_temp2(ntr),0,[-1000,3500]);
                    t2 = tlo2:thi2;
                    rate2 = 1000*smooth_es(spks2', wid, noedge)';
                    rate_opp_cue(ntr,:) = rate2;
                end
                for nt = 1:size(rate_cue,2)
                    neuron_ROC(n,nt) = arrayROC(rate_cue(:,nt),rate_opp_cue(:,nt));
                end
            catch
                disp(['error processing neuron  ', num2str(Neurons(n)) '  Dir1=' num2str(Best_class(n))])
            end
        end
        NN1 = size(neuron_ROC,1);
        mr1 = mean(neuron_ROC, 1); % mean rate
        semr1 = std(neuron_ROC, 1)/sqrt(NN1); % SEM
        plt_save.t{g,s} = t1;
        plt_save.Nt{g,s} = NN1;
        plt_save.mr{g,s} = mr1;
        plt_save.semr{g,s} = semr1;
        plt_save.roc{g,s} = neuron_ROC;
    end
end
%% time resolved permutation test
p_value = [];
for g = 1:size(group,1)
    for nt = 1:size(t1,2)
        % Prepare data
        data1 = plt_save.roc{g,1}(:,nt);
        data2 = plt_save.roc{g,2}(:,nt);
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
        if mod(nt,100)==0
            disp([num2str(nt)])
        end
        % Compute the p-value
        p_value(g,nt) = sum(abs(perm_diffs) >= abs(obs_diff)) / num_permutations;
    end
end