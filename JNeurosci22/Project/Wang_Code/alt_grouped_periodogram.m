f_for_plot =( lb:f_interval:ub) + f_interval/2;
fig_size = [5 1 8 5];
fig_lib_temp = 'fig_temp';
%%  yp-ya
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram_std(periodos_young_ps, f_for_plot, [.0, .0, .9]);
lines(2) = zw_plot_periodogram_std(periodos_young_as, f_for_plot, [.1, .1, .4]);
legend(lines, {'Young-PS', 'Young-AS'});
ylim([0, 0.14]);
print(gcf, fullfile(project_dir, fig_lib_temp, 'yp-ya'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib_temp, 'yp-ya'), '-dtiff', '-r400')
%%  ap-aa
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram_std(periodos_adult_ps, f_for_plot, [.9, .0, .0]);
lines(2) = zw_plot_periodogram_std(periodos_adult_as, f_for_plot, [.4, .1, .1]);
legend(lines, {'Adult-PS', 'Adult-AS'});
ylim([0, 0.14]);
print(gcf, fullfile(project_dir, fig_lib_temp, 'ap-aa'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib_temp, 'ap-aa'), '-dtiff', '-r400')
%%  all
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram_std(periodos_adult_ps, f_for_plot, [.9, .0, .0]);
lines(2) = zw_plot_periodogram_std(periodos_adult_as, f_for_plot, [.4, .1, .1]);
lines(3) = zw_plot_periodogram_std(periodos_young_ps, f_for_plot, [.0, .0, .9]);
lines(4) = zw_plot_periodogram_std(periodos_young_as, f_for_plot, [.1, .1, .4]);
significant = f_for_plot(f1(:, 2) < .05);
for i = 1:length(significant)
    lines(5) = plot(significant(i) + [-f_interval/2 ,0, f_interval/2], zeros(1, 3) + 0.01, 'Color', [.6, .6, .6], 'LineWidth', 4);
end
legend(lines, {'Adult-PS', 'Adult-AS', 'Young-PS', 'Young-AS', 'Significant stage main effect'});
ylim([0, 0.14]);
xlim([0, 60])
print(gcf, fullfile(project_dir, fig_lib_temp, 'all'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib_temp, 'all'), '-dtiff', '-r400')
%%  ap-yp
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram_std(periodos_adult_ps, f_for_plot, [.9, .0, .0]);
lines(2) = zw_plot_periodogram_std(periodos_young_ps, f_for_plot, [.0, .0, .9]);
legend(lines, {'Adult-PS', 'Young-PS'});
ylim([0, 0.14]);
print(gcf, fullfile(project_dir, fig_lib_temp, 'ap-yp'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib_temp, 'ap-yp'), '-dtiff', '-r400')
%%  aa-ya
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram_std(periodos_adult_as, f_for_plot, [.4, .1, .1]);
lines(2) = zw_plot_periodogram_std(periodos_young_as, f_for_plot, [.1, .1, .4]);
legend(lines, {'Adult-AS', 'Young-AS'});
ylim([0, 0.14]);
print(gcf, fullfile(project_dir, fig_lib_temp, 'aa-ya'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib_temp, 'aa-ya'), '-dtiff', '-r400')
%%