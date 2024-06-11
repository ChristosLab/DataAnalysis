%   Pre-allocation
title_strings = {'All sites', 'Young sites', 'Adult sites'};
gps = {[site_tuning_cat{1:2, [1,3,4], :, 2}], [site_tuning_cat{1:2, [1,3,4], :, 1}]; ...
    [site_tuning_cat{1, [1,3,4], :, 2}], [site_tuning_cat{1, [1,3,4], :, 1}]; ...
    [site_tuning_cat{2, [1,3,4], :, 2}], [site_tuning_cat{2, [1,3,4], :, 1}]};
stat = cell(0);

%%  Bootstrapping inform. vs non-inform.
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
fname_ = 'lfp_pev_cluster_stat';
save(fullfile(project_dir, output_database, fname_), 'stat');