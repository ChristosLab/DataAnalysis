% neuron firing rate analysis of PFC neurons at diff layers in ODR task
% For ODR task, sig neurons, baseline rate
% 20230619, Junda Zhu
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
%% compute baseline
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
for g = 1:size(group,1)
    nn = [];
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_Cue = neuron_info.cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        fix_rate = [];
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
            fix_rate(n) = mean(rate_fix,"omitnan");
        end
        fix_rate_mean = mean(fix_rate,'omitnan');
    fix_rate_sem = std(fix_rate,'omitnan')./sqrt(length(fix_rate)); 
    nn = size(odr_data_group,1);
        % save for plot
        plt_save.fix_rate{g,s} = fix_rate;
        plt_save.fix_rate_sem{g,s} = fix_rate_sem;
        plt_save.nn{g,s} = nn;
        plt_save.fix_rate_mean{g,s} = fix_rate_mean;
    end
end
%% t test
for g = 1:size(group,1)
    [h(g),p(g),ci,stat] = ttest2(plt_save.fix_rate{g,1},plt_save.fix_rate{g,2})
end