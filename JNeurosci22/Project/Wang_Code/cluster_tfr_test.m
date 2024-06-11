%%
cfg = [];
% cfg.elec             = elec;
% cfg.neighbourdist    = 4;
cfg.latency          = 'all';
cfg.frequency        = 'all';
% cfg.channel          = 'EEG1010' % see FT_CHANNELSELECTION
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
cfg.avgoverchan      = 'no';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.numrandomization = 200;
cfg.correctm         = 'cluster';
cfg.method           = 'montecarlo';
cfg.design           = [
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10     % subject number
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2  ];  % condition number
cfg.uvar = 1;                                   % "subject" is unit of observation
cfg.ivar = 2;                                   % "condition" is the independent variable
ld = [];
ld.powspctrm = randn([10, 1, 100,100]);
ld.dimord = 'chan_freq_time';
ld.freq   = 1:100;
ld.label = {'1'};
ld.time   = (1:100)./1000;
sd = ld;
sd.powspctrm = randn([10, 1, 100,100]);
stat = ft_freqstatistics(cfg, ld, sd)
%%
nb_cfg = [];
nb_cfg.method = 'distance';
nb_cfg.channel = 'all';
nb = ft_prepare_neighbours(nb_cfg, ld)
%%
cfg = [];
cfg.alpha = 0.025;
elec = [];
elec.chanpos = [0 0 0];
elec.label = {'1'};
stat.elec = elec;
ft_clusterplot(cfg, stat)
%%
x_range = 1:1:20;
y_range = 1:1:20;
test_mat = eye(numel(y_range), numel(x_range));
test_mat = zeros(numel(y_range), numel(x_range));
test_mat(10:15,5:8) = 2;
test_mat(6:8,3:4) = 1;
test_mat(8, 4) = 0;
[test_x, test_y] = meshgrid(x_range, y_range);
figure;
hold on
% subplot(1,2,1);
imagesc(x_range, y_range, test_mat);
% subplot(1,2,2);
[test_c, test_h] = contour(test_x, test_y, test_mat, 1);
test_cc = contourc(x_range, y_range, test_mat, [.5, .5]);
test_ccc = zw_plot_cluster_contour(test_mat);
for i = 1:numel(test_ccc)
    plot(test_ccc{i}(1, :), test_ccc{i}(2, :), '-r', 'LineWidth', 3)
end
%%
cfg = [];
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
cfg.avgoverchan      = 'no';
cfg.statistic        = 'ft_statfun_indepsamplesT';
% cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.numrandomization = 10000;
cfg.tail             = 0;
cfg.correctm         = 'cluster';
cfg.method           = 'montecarlo';
% 
gp_y = [site_mod_cat{1, : , 1:2}];
gp_a = [site_mod_cat{2, : , 1:2}];
gp_y_i = [site_tuning_cat{1, [1,3,4], :, 2}];
gp_a_i = [site_tuning_cat{2, [1,3,4], :, 2}];
gp_y_n = [site_tuning_cat{1, [1,3,4], :, 1}];
gp_a_n = [site_tuning_cat{2, [1,3,4], :, 1}];

%
in_y = [];
in_y.dimord = 'chan_freq_time';
in_y.freq   = 2:2:128;
% in_y.freq   = 2:2:64;
in_y.label = {'1'};
% in_y.time   = [1/50:1/50:5] - 1;
% in_y.time   = [1/50:1/50:4];
in_y.time   = [1/50:1/50:2];
%
in_a = in_y;
% %
% powspctrm_y = temp_cwt_cue_b(gp_y, :, :);
% powspctrm_a = temp_cwt_cue_b(gp_a, :, :);
powspctrm_y = temp_cwt_cue_b(gp_y, :, 51:150);
powspctrm_a = temp_cwt_cue_b(gp_a, :, 51:150);
% powspctrm_y = temp_cwt_cue_b(gp_y_i, :, 51:150);
% powspctrm_a = temp_cwt_cue_b(gp_a_i, :, 51:150);
% powspctrm_y = temp_cwt_cue_b(gp_y_n, :, 51:150);
% powspctrm_a = temp_cwt_cue_b(gp_a_n, :, 51:150);
% powspctrm_y = temp_cwt_cue_b(gp_y, :, 51:end);
% powspctrm_a = temp_cwt_cue_b(gp_a, :, 51:end);
% powspctrm_y = temp_cwt_cue_b(gp_y, 1:32, :);
% powspctrm_a = temp_cwt_cue_b(gp_a, 1:32, :);
% powspctrm_y = temp_cwt_cue_b(gp_y, 1:32, 51:end);
% powspctrm_a = temp_cwt_cue_b(gp_a, 1:32, 51:end);
powspctrm_y = reshape(powspctrm_y, [size(powspctrm_y, 1), 1, size(powspctrm_y, 2), size(powspctrm_y, 3)]);
powspctrm_a = reshape(powspctrm_a, [size(powspctrm_a, 1), 1, size(powspctrm_a, 2), size(powspctrm_a, 3)]);
%
in_y.powspctrm = powspctrm_y;
in_a.powspctrm = powspctrm_a;
%
cfg.design = [1:(numel(gp_y) + numel(gp_a)); ones(1, numel(gp_y)), 1 + ones(1, numel(gp_a))];
% cfg.design = [1:(numel(gp_y_i) + numel(gp_a_i)); ones(1, numel(gp_y_i)), 1 + ones(1, numel(gp_a_i))];
% cfg.design = [1:(numel(gp_y_n) + numel(gp_a_n)); ones(1, numel(gp_y_n)), 1 + ones(1, numel(gp_a_n))];
% cfg.uvar = 1;                                   % "subject" is unit of observation
cfg.ivar = 2;                                   % "condition" is the independent variable
% cfg.design           = [
%   1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10     % session X electrode (trivial observations)
%   1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2  ];  % Developmental stage
stat = ft_freqstatistics(cfg, in_y, in_a);
%%
figure;
hold on
imagesc(squeeze(stat.prob));
%%
out_vertices = zw_plot_cluster_contour(squeeze(stat.posclusterslabelmat));
%%
figure;
hold on
imagesc(squeeze(stat.prob));
out_vertices = zw_plot_cluster_contour(squeeze(stat.negclusterslabelmat));
for i = 1:numel(out_vertices)
    plot(out_vertices{i}(1, :), out_vertices{i}(2, :), '-r', 'LineWidth', 3)
end
%%
for_plot_pos = squeeze(stat.posclusterslabelmat);
for_plot_neg = squeeze(stat.negclusterslabelmat);
for_plot     = zeros(size(for_plot_pos));
for i = 1:numel(for_plot_pos)
    if for_plot_pos(i)
        for_plot(i) = stat.posclusters(for_plot_pos(i)).prob <0.05;
    end
    if for_plot_neg(i)
        for_plot(i) = -(stat.negclusters(for_plot_neg(i)).prob <0.05);
    end
end
figure('Units', 'inches', 'Position', [2, 2, 4, 3.2]);
hold on
imagesc(in_y.time, in_y.freq, for_plot, [-1,1]);
line()
xlabel('Time from stim. onset (s)')
ylabel('Frequency (Hz)')
xlim([in_y.time(1), in_y.time(end)]);
ylim([in_y.freq(1), in_y.freq(end)]);
xticks([0:0.5:2])
set(gca, {'FontSize', 'FontWeight'}, {11, 'bold'})
colormap([0, 0, .8; 0, 0, 0; .8, 0, 0]);
h = colorbar;
h.Ticks = [-.7, 0, .7];
h.TickLabels = {'Adolescent\newline<\newlineAdult', 'Non-sig.', 'Adolescent\newline>\newlineAdult'};
h.FontSize = 7;
%%
figure;
hold on
imagesc(squeeze(stat.mask));
contour(squeeze(stat.mask), 1);
%%
figure;
hold on
imagesc(squeeze(stat.negclusterslabelmat));
%%
figure;
hold on
imagesc(squeeze(stat.posclusterslabelmat));
%%
stat_all = {stat, cfg};
%%
stat_young = stat;
%%
figure; hold on
plot(2:2:128, mean(mean(temp_cwt_cue(gp_y,:,1:50), 3), 1))
plot(2:2:128, mean(mean(temp_cwt_cue(gp_a,:,1:50), 3), 1))
%%
figure; hold on
plot(2:2:128, mean(mean(temp_cwt_cue_b(gp_y,:,76:150), 3), 1), 'r')
plot(2:2:128, mean(mean(temp_cwt_cue_b(gp_a,:,76:150), 3), 1), 'b')

plot(2:2:128, mean(mean(temp_cwt_cue_b(gp_y,:,51:75), 3), 1), '--r')
plot(2:2:128, mean(mean(temp_cwt_cue_b(gp_a,:,51:75), 3), 1), '--b')

%%
figure; hold on
plot(2:2:128, mean(log(mean(temp_cwt_cue(gp_y,:,26:50), 3)), 1))
% plot(2:2:128, mean(mean(temp_cwt_cue(gp_y,:,51:75), 3), 1))
plot(2:2:128, mean(log(mean(temp_cwt_cue(gp_y,:,76:150), 3)), 1))

%%
figure; hold on
plot(2:2:128, mean(log(mean(temp_cwt_cue(gp_a,:,26:50), 3)), 1))
% plot(2:2:128, mean(mean(temp_cwt_cue(gp_y,:,51:75), 3), 1))
plot(2:2:128, mean(log(mean(temp_cwt_cue(gp_a,:,76:150), 3)), 1))
%%
function find_cluster_center(labelmat)

end