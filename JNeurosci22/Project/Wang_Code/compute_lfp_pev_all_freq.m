%%  lfp class ANOVA
%%
zw_setpath;
%%
fname_ = 'complete_cwt_repo_3_2_2-21';
load(fullfile(project_dir, output_database, fname_));
%%
target_frs = [4,8; 8,16; 16, 32; 32, 64];
n_bands = size(target_frs, 1);
baseline_bin = [26, 50];
%%  CWT parameters
fs = 500;
cue_dur = [-1 + 1/fs, 4]; % Signal length 
normalizer = 4;
kernel_flag = 3;
down_sample = 10;
f_range = 2:2:128;
avg_method = 1; %   0 - arithmetic; 1 - geometric
%%
baseline_normalized_tfr_repo_all_freq = struct();
baseline_normalized_tfr_repo_all_freq.class.cue_tfr = cell(numel(f_range), 1);
% baseline_normalized_tfr_repo_all_freq.class.sac_tfr = cell(8,1);
%%  Compute tfr
for i = [t{:}]
    tic
    n_class_ = numel(cwt_repo(i).class);
    cue_cwt_ = {cwt_repo(i).class.cue_cwt};
%     sac_cwt_ = {cwt_repo(i).class.sac_cwt};
    for k = 1:numel(f_range)
            for m = 1:n_class_
               temp_tfr_ = sum(cwt_repo(i).class(m).cue_cwt(:, k, :), 2)...
                   ./nanmean(sum(cwt_repo(i).class(m).cue_cwt(:, k, baseline_bin(1):baseline_bin(2)), 2), 3);
               temp_tfr_ = reshape(temp_tfr_, [numel(temp_tfr_)/250, 250]);
                baseline_normalized_tfr_repo_all_freq(i).class(m).cue_tfr{k} = temp_tfr_;
%                 baseline_normalized_tfr_repo_all_freq(i).class(m).sac_tfr{k} = squeeze(sum(cwt_repo(i).class(m).sac_cwt(:, target_frs(k, 1):target_frs(k, 2), :), 2));
            end
            %             [~, out_tb_, out_stats_] = anova1(anova_y_, anova_label_, 'off');
            %             neuron_analysis_output_cue(i, j).tb = out_tb_;
            %             neuron_analysis_output_cue(i, j).stats = out_stats_;
    end
    toc
end
%%  ANOVA and PEV
n_bins_lfp = (diff(cue_dur)*fs + 1)/down_sample;
lfp_pev_cue_all_freq = [];
for i = [t{:}]
    tic
    n_class_ = numel(baseline_normalized_tfr_repo_all_freq(i).class);
    cue_tfr_ = {baseline_normalized_tfr_repo_all_freq(i).class.cue_tfr};
%     sac_cwt_ = {cwt_repo(i).class.sac_cwt};
    for k = 1:numel(f_range)
        for j = 1:n_bins_lfp
            anova_y_cue_ = [];
            anova_label_cue_ = [];
%             anova_y_sac_ = [];
%             anova_label_sac_ = [];
            for m = 1:n_class_
                if m > 8
                    class_id_ = mod(m, 8);
                else
                    class_id_ = m;
                end
                anova_y_cue_ = [anova_y_cue_, cue_tfr_{m}{k}(:, j)'];
                anova_label_cue_ = [anova_label_cue_, class_id_ + zeros(size(cue_tfr_{m}{k}(:, j)'))];
            end
            [~, out_tb_, ~] = anova1(anova_y_cue_, anova_label_cue_, 'off');
            try
                lfp_pev_cue_all_freq(i, j, k) = zw_pev(out_tb_);
            catch
                i
                lfp_pev_cue_all_freq(i, j, k) = nan;
            end

        end
    end
    toc
end
%%
lfp_pev_cue_all_freq = permute(lfp_pev_cue_all_freq, [1, 3, 2]);
%%
fname_ = 'lfp_pve_baseline_normalized_dur_all_freq.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_pev_cue_all_freq');
%%  Checking plots
%
figure();
imagesc([1/50-1:1/50:4], f_range, squeeze(mean(lfp_pev_cue_all_freq([site_tuning_cat{1, :, 1, :}], :, :), 1)), [0, 0.1])
set(gca, 'YDir', 'normal');
title('Adolescent non-informative');
%
figure();
imagesc([1/50-1:1/50:4], f_range, squeeze(mean(lfp_pev_cue_all_freq([site_tuning_cat{1, :, 2, :}], :, :), 1)), [0, 0.1])
set(gca, 'YDir', 'normal');
title('Adolescent informative');
%
figure();
imagesc([1/50-1:1/50:4], f_range, squeeze(mean(lfp_pev_cue_all_freq([site_tuning_cat{2, :, 1, :}], :, :), 1)), [0, 0.1])
set(gca, 'YDir', 'normal');
title('Adult non-informative');
%
figure();
imagesc([1/50-1:1/50:4], f_range, squeeze(mean(lfp_pev_cue_all_freq([site_tuning_cat{2, :, 2, :}], :, :), 1)), [0, 0.1])
set(gca, 'YDir', 'normal');
title('Adult informative');
%%  Cluster-based analysis
%   Fieldtrip path-setting
addpath(fullfile(project_dir, code_lib, 'External\fieldtrip\'));
ft_defaults
%   Pre-allocation
title_strings = {'All sites', 'Informative sites', 'Non-informative sites'};
gps = {[site_mod_cat{1, : , 1:2}], [site_mod_cat{2, : , 1:2}]; ...
    [site_tuning_cat{1, [1,3,4], :, 2}], [site_tuning_cat{2, [1,3,4], :, 2}]; ...
    [site_tuning_cat{1, [1,3,4], :, 1}], [site_tuning_cat{2, [1,3,4], :, 1}]};
stat = cell(0);

%%  Bootstrapping y. vs a.
for i = 1:numel(title_strings)
    gp_y = gps{i, 1};
    gp_a = gps{i, 2};
    cfg = [];
    cfg.latency          = 'all';
    cfg.frequency        = 'all';
    cfg.avgovertime      = 'no';
    cfg.avgoverfreq      = 'no';
    cfg.avgoverchan      = 'no';
    cfg.statistic        = 'ft_statfun_indepsamplesT';
    cfg.numrandomization = 10000;
    cfg.tail             = 0;
    cfg.correctm         = 'cluster';
    cfg.method           = 'montecarlo';
    %
    %
    in_y = [];
    in_y.dimord = 'chan_freq_time';
    in_y.freq   = 2:2:128;
    in_y.label = {'1'};
    in_y.time   = [1/50:1/50:2];
    in_y.time   = [1/50:1/50:3];
    in_a = in_y;
%     pevspctrm_y = lfp_pev_cue_all_freq(gp_y, :, 51:150);
%     pevspctrm_a = lfp_pev_cue_all_freq(gp_a, :, 51:150);
    pevspctrm_y = lfp_pev_cue_all_freq(gp_y, :, 51:200);
    pevspctrm_a = lfp_pev_cue_all_freq(gp_a, :, 51:200);
    pevspctrm_y = reshape(pevspctrm_y, [size(pevspctrm_y, 1), 1, size(pevspctrm_y, 2), size(pevspctrm_y, 3)]);
    pevspctrm_a = reshape(pevspctrm_a, [size(pevspctrm_a, 1), 1, size(pevspctrm_a, 2), size(pevspctrm_a, 3)]);
    %
    in_y.powspctrm = pevspctrm_y;
    in_a.powspctrm = pevspctrm_a;
    %
    cfg.design = [1:(numel(gp_y) + numel(gp_a)); ones(1, numel(gp_y)), 1 + ones(1, numel(gp_a))];
    cfg.ivar = 2;                                   % "condition" is the independent variable
    stat{i} = ft_freqstatistics(cfg, in_y, in_a);
end
%%
%%  Plotting
for i = 1:numel(title_strings)
    for_plot_pos = squeeze(stat{i}.posclusterslabelmat);
    for_plot_neg = squeeze(stat{i}.negclusterslabelmat);
    for_plot     = zeros(size(for_plot_pos));
    for j = 1:numel(for_plot_pos)
        if for_plot_pos(j)
            for_plot(j) = stat{i}.posclusters(for_plot_pos(j)).prob <0.05;
        end
        if for_plot_neg(j)
            for_plot(j) = -(stat{i}.negclusters(for_plot_neg(j)).prob <0.05);
        end
    end
    figure('Units', 'inches', 'Position', [2, 2, 4, 3.2]);
    hold on
    imagesc(in_y.time, in_y.freq, for_plot, [-1,1]);
    pos_center = zw_find_cluster_center(for_plot == 1);
    if ~any(isnan(pos_center))
        text(in_y.time(pos_center(2)), in_y.freq(pos_center(1)), sprintf('%.3f', stat{i}.posclusters(1).prob), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
    neg_center = zw_find_cluster_center(for_plot == -1);
    if ~any(isnan(neg_center))
        text(in_y.time(neg_center(2)), in_y.freq(neg_center(1)), sprintf('%.3f', stat{i}.negclusters(1).prob), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
    plot([in_y.time(1), in_y.time(end)], [1, 1].*32, '--', 'Color', 'w')
    xlabel('Time from stim. onset (s)')
    ylabel('Frequency (Hz)')
    xlim([in_y.time(1), in_y.time(end)]);
    ylim([in_y.freq(1), in_y.freq(end)]);
    xticks([0:0.5:2])
    yt = yticks;
    yticks(sort([yt, 32]));
    title(title_strings{i})
    set(gca, {'FontSize', 'FontWeight'}, {11, 'bold'})
    colormap([0, 0, .8; 0, 0, 0; .8, 0, 0]);
    h = colorbar;
    h.Ticks = [-.7, 0, .7];
    h.TickLabels = {'Adolescent\newline<\newlineAdult', 'Non-sig.', 'Adolescent\newline>\newlineAdult'};
    h.FontSize = 7;
end
%%
stat_0_to_2 = stat;
%%
stat_0_to_3 = stat;
%%
