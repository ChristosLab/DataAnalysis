% neuron tuning analysis of PFC neurons at diff layers in ODR task
% For ODR task, sig neurons, width (std of the gaussian curve)
% generate .csv files for ANOVA using python (under the same folder)
% 20230710, Junda Zhu
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
%% compute firing rate
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
for g = 1:size(group,1)
    nn = [];
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_class = neuron_info.cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        rate_test = [];
        rate_test_sort = [];
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
            rate_test(n,:) = rate_del;
            [~, max_test_class] = max(rate_test(n,:));
            % align max firing rate
            for n_class = 1:8
                if n_class >= max_test_class
                    rate_test_sort(n,n_class-max_test_class+1) = rate_test(n,n_class);
                else
                    rate_test_sort(n,n_class-max_test_class+1+8) = rate_test(n,n_class);
                end
            end
        end
        rate_test_mean = mean(rate_test_sort,'omitnan');
        rate_test_sem = std(rate_test_sort,'omitnan')./sqrt(length(rate_test_sort));
        nn = size(odr_data_group,1);
        % save for plot
        plt_save.rate{g,s} = rate_test_sort;
        plt_save.rate_sem{g,s} = rate_test_sem;
        plt_save.nn{g,s} = nn;
        plt_save.rate_mean{g,s} = rate_test_mean;
    end
end
%% construct table
tbl_for_anova = table;
for ii = 1:size(plt_save.rate,1)
    for jj = 1:size(plt_save.rate,2)
        tbl_temp = table;
        tbl_temp.rate = plt_save.rate{ii,jj}(:);
        loc = repmat([1,2,3,4,5,6,7,8],size(plt_save.rate{ii,jj},1),1);
        tbl_temp.loc = loc(:);
        tbl_temp.group = ones(numel(plt_save.rate{ii,jj}),1)*ii;
        tbl_temp.stage = ones(numel(plt_save.rate{ii,jj}),1)*jj;
        tbl_for_anova = [tbl_for_anova;tbl_temp];
    end
end
writetable(tbl_for_anova,'tbl_for_anova_del.csv');