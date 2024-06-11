%%
clm_n = 100;
clm_ = [[zeros(1, clm_n);zeros(1, clm_n); linspace(1, 0, clm_n)]'; [linspace(0, 1, clm_n); zeros(1, clm_n);zeros(1, clm_n)]'];

% clm_ = [[zeros(1, clm_n);zeros(1, clm_n); log10(linspace(10, 1, clm_n))]'; [log10(linspace(1, 10, clm_n)); zeros(1, clm_n);zeros(1, clm_n)]'];
%%  Band power change
%% 
figure(); plot(mean(squeeze(nanmean(temp_cwt(t14, 12:17, :)./temp_baseline(t14, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t24, 12:17, :)./temp_baseline(t24, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt(t14, 32:64, :)./temp_baseline(t14, 32:64), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t24, 32:64, :)./temp_baseline(t24, 32:64), 1)), 1));%ylim([6, 8]*10^-3)
figure(); plot(mean(squeeze(nanmean(temp_cwt(t14, 8:16, :)./temp_baseline(t14, 8:16), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t24, 8:16, :)./temp_baseline(t24, 8:16), 1)), 1));%ylim([6, 8]*10^-3)
%%  Modulation depth distribution
%%
hist(nanmean(nanmean(temp_cwt(t24, 8:16, 51:75)./temp_baseline(t24, 8:16), 2), 3))
%%
hist(nanmean(nanmean(temp_cwt(t24,32:64, 51:75)./temp_baseline(t24, 32:64), 2), 3), 0.5:0.1:2)
%%
figure()
hist(nanmean(nanmean(temp_cwt(t24,16:32, 51:75)./temp_baseline(t24, 16:32), 2), 3), 0.5:0.1:2)
%%
figure()
hist(nanmean(nanmean(temp_cwt(t13,25:64, 51:75)./temp_baseline(t13, 25:64), 2), 3), 0.5:0.1:2)
%%
figure()
histogram(nanmean(nanmean(temp_cwt(t14, 4:8, 51:75)./temp_baseline(t14, 4:8), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.1:2)
hold on

histogram(nanmean(nanmean(temp_cwt(t24, 4:8, 51:75)./temp_baseline(t24, 4:8), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.1:2)
%%
figure()
histogram(nanmean(nanmean(temp_cwt(t14, 8:16, 51:75)./temp_baseline(t14, 8:16), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
hold on

histogram(nanmean(nanmean(temp_cwt(t24, 8:16, 51:75)./temp_baseline(t24, 8:16), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
legend({'Adolescent', 'Adult'})
%%
figure()
histogram(nanmean(nanmean(temp_cwt(t14, 32:64, 51:75)./temp_baseline(t14, 32:64), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
hold on

histogram(nanmean(nanmean(temp_cwt(t24, 32:64, 51:75)./temp_baseline(t24, 32:64), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
legend({'Adolescent', 'Adult'})
%%
figure();
imagesc(squeeze(nanmean(temp_cwt_b(t14, :, :), 1)) - 1, [-.4,.4]);
set(gca, 'YDir','normal'); 
colormap(clm_)
h = colorbar();
%%
figure();
imagesc([-1, 3], [2,128], (squeeze(nanmean(temp_cwt_b([t13;t14], :, :), 1)) - 1)*100, [-.4,.4]*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, 'Power (% change from fixation)')
set_plot_12_4();
title('Adolescent')
fig_name = 'spectrogram_change_young'
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
% 
figure();
imagesc([-1, 3], [2,128], (squeeze(nanmean(temp_cwt_b([t23;t24], :, :), 1)) - 1)*100, [-.4,.4]*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, 'Power (% change from fixation)')
set_plot_12_4();
title('Adult')
fig_name = 'spectrogram_change_adult'
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');