% neuron tuning analysis of PFC neurons at diff layers in ODR task
% For ODR task, sig neurons, width (std of the gaussian curve)
% calculate firing rate in epochs from correct inRF trials pooled
% plot grouped tuning curves
% 20230619, Junda Zhu
%% load data
clearvars
load('sig_odr_data_depth_20230906_raw_max.mat');
file_name = string(neuron_info.Filename);
neuron_info.ID = extract(file_name,1);
%% seg data/label groups
neuron_info.group(neuron_info.Depth<=800)=1;
neuron_info.group(neuron_info.Depth>800&neuron_info.Depth<=1200)=2;
neuron_info.group(neuron_info.Depth>1200)=3;
odr_data = odr_data(~isnan(neuron_info.Depth),:);
neuron_info = neuron_info(~isnan(neuron_info.Depth),:);
%% select neurons
odr_data = odr_data(neuron_info.sulcus==0,:);
neuron_info = neuron_info(neuron_info.sulcus==0,:);
%% compute tuning curve
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
for g = 1:size(group,1)
    nn = [];
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_Cue = neuron_info.del(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        class_mean_all = [];
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
            rate_test = rate_cue; % which to look at
            [~, max_test_class] = max(rate_test(1:8));
            %[a max_class] = min(varCue(1:8));   % for inhibit neurons

            % align max firing rate
            for n_class = 1:8
                if n_class >= max_test_class
                    class_mean(n_class-max_test_class+1) = rate_test(n_class);
                else
                    class_mean(n_class-max_test_class+1+8) = rate_test(n_class);
                end
            end
            class_mean_order = [5,4,3,2,1,8,7,6];
            class_mean_all(1:8,n) = class_mean(class_mean_order);
            class_mean_all(9,n) = class_mean_all(1,n);
        end
        % gaussian fit
        [st, m, sem, loctt, rateyy, d, R] = gaus_fit_8_loc(class_mean_all);
        fix_rate_mean_all = mean(fix_rate,'omitnan');
        nn = size(odr_data_group,1);
        % save for plot
        plt_save.st{g,s} = st;
        plt_save.m{g,s} = m;
        plt_save.sem{g,s} = sem;
        plt_save.d{g,s} = d;
        plt_save.loctt{g,s} = loctt;
        plt_save.rateyy{g,s} = rateyy;
        plt_save.nn{g,s} = nn;
        plt_save.fix_rate_mean_all{g,s} = fix_rate_mean_all;
    end
end
%% plot fitting curve
% define a few of colors
my_color = linspecer(3);
my_color = my_color([1:3],:);

figure
clf
set(gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.2,0.1,0.6,0.8] );
for g = 1:size(group,1)
    sp = subplot(1,3,g);
    for s = 1:size(stage,1)
        hold on
        % plot curves
        plot(plt_save.loctt{g,s}, plt_save.rateyy{g,s}, 'Color',my_color(g,:)/s, 'linewidth', 4)
        % plot errors as a shaded area
        errorbar(plt_save.st{g,s},plt_save.m{g,s},plt_save.sem{g,s},'marker','o','linestyle','none','Color',my_color(g,:)/s,HandleVisibility='off')
        fix_edges=1:.5:9;
        line(fix_edges,plt_save.fix_rate_mean_all{g,s}'*ones(1,length(fix_edges)),'Color',my_color(g,:)/s,'linestyle','--','linewidth', 2,HandleVisibility='off');
    end
    axis([0.5 9.5 0 30])
    xlabel('Spatial location')
    ylabel('firing rate (sp/s)')
    title('tuning curve')
    legend({'adult','young'})
    annotation('textbox','Position',sp.Position,'Vert','bottom','FitBoxToText','on','String', "std tuning curve: "+[plt_save.d{g,:}]+" n = "+[plt_save.nn{g,:}])
end
sgtitle('delay rate');
%% additional: compute firing rate from specific trial epoch
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
epoch_start = 0.5;
epoch_end= 2.0;
for g = 1:size(group,1)
    nn = [];
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_Cue = neuron_info.cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        for n = 1:length(Neurons)
            MatData = odr_data_group(n,:);
            rate_test = [];
            for j = 1:length(MatData)
                try
                    [FR_temp, ntrs_temp] = Get_FRbyneuron_AllTrials_alignCue(MatData,j,epoch_start,epoch_end);
                    rate_test(j) = mean(FR_temp);
                    if isfield(MatData{j},'fixrate')
                        rate_fix(j) = mean([MatData{j}.fixrate]);
                    else
                        rate_fix(j) = mean([MatData{j}.fix]);
                    end
                catch
                    rate_test(j) = nan;
                end
            end
            fix_rate(n) = mean(rate_fix,"omitnan");
            [~, max_test_class] = max(rate_test(1:8));
            %[a max_class] = min(varCue(1:8));   % for inhibit neurons

            % align max firing rate
            for n_class = 1:8
                if n_class >= max_test_class
                    class_mean(n_class-max_test_class+1) = rate_test(n_class);
                else
                    class_mean(n_class-max_test_class+1+8) = rate_test(n_class);
                end
            end
            class_mean_order = [5,4,3,2,1,8,7,6];
            class_mean_all(1:8,n) = class_mean(class_mean_order);
            class_mean_all(9,n) = class_mean_all(1,n);
        end
        % gaussian fit
        [st, m, sem, loctt, rateyy, d, R] = gaus_fit_8_loc(class_mean_all);
        fix_rate_mean_all = mean(fix_rate,'omitnan');
        nn = size(odr_data_group,1);
        % save for plot
        plt_save.st{g,s} = st;
        plt_save.m{g,s} = m;
        plt_save.sem{g,s} = sem;
        plt_save.d{g,s} = d;
        plt_save.loctt{g,s} = loctt;
        plt_save.rateyy{g,s} = rateyy;
        plt_save.nn{g,s} = nn;
        plt_save.fix_rate_mean_all{g,s} = fix_rate_mean_all;
    end
end