%%
last_young_session = [20, 26, 66, 27];
%   Convert  dummy code
[~, last_young_session] = max(last_young_session == sessions', [], 1);
%%
monk_ind = 3;
spectrogram_change_young_ps = [];
spectrogram_change_young_as = [];
spectrogram_change_adult_ps = [];
spectrogram_change_adult_as = [];
tic
for a = 1:size(groups, 1)
    for b = 1:size(groups, 2)
        for c = 1:size(groups, 3)
            c
%             for d = 1:size(groups, 4)
            for d = monk_ind
                p_ = squeeze(spectrogram_change(a, b, c, d, :, :));
                if any(isnan(p_))
                    continue
                end
                if c > last_young_session(d)
                    if b == 1
                        spectrogram_change_adult_ps = cat(3, spectrogram_change_adult_ps, p_);
                    else
                        spectrogram_change_adult_as = cat(3, spectrogram_change_adult_as, p_);
                    end
                else
                    if b == 1
                        spectrogram_change_young_ps = cat(3, spectrogram_change_young_ps, p_);
                    else
                        spectrogram_change_young_as = cat(3, spectrogram_change_young_as, p_);
                    end
                end
            end
            toc
        end
    end
end
%%
temp_input = cell(4, 1);
temp_input{1} = permute(spectrogram_change_young_ps, [3, 1, 2]); 
temp_input{2} = permute(spectrogram_change_young_as, [3, 1, 2]); 
temp_input{3} = permute(spectrogram_change_adult_ps, [3, 1, 2]); 
temp_input{4} = permute(spectrogram_change_adult_as, [3, 1, 2]);
%%
tic
% tbl_temp = [];
ttest_temp = [];
for ind = 1:(size(temp_input{1}, 2)*size(temp_input{1}, 3))
%     tbl_temp = cat(3, tbl_temp, zw_quick_anova(temp_input{1}, temp_input{2}, temp_input{3}, temp_input{4}, ind));
    ttest_temp = cat(1, ttest_temp, zw_ttest2(temp_input{1}, temp_input{2}, temp_input{3}, temp_input{4}, ind));
    if mod(ind, 100) == 0
        ind
        toc
    end
end
%%  Save group output by monkey
fname_ = sprintf('spectrogram_change_stat_%d_%d_%d_%d_monkey_%d).mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2), monk_ind);
save(fullfile(project_dir, output_database, fname_), 'tbl_temp', 'monk_ind');

fname_ = sprintf('spectrogram_change_input_%d_%d_%d_%d_monkey_%d).mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2), monk_ind);
save(fullfile(project_dir, output_database, fname_), 'temp_input', 'monk_ind');

fname_ = sprintf('spectrogram_change_t_%d_%d_%d_%d_monkey_%d).mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2), monk_ind);
save(fullfile(project_dir, output_database, fname_), 'ttest_temp', 'monk_ind');

%%  Load saved group output by monkey
monk_ind = 3;
fname_ = sprintf('spectrogram_change_stat_%d_%d_%d_%d_monkey_%d).mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2), monk_ind);
load(fullfile(project_dir, output_database, fname_), 'tbl_temp');

fname_ = sprintf('spectrogram_change_input_%d_%d_%d_%d_monkey_%d).mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2), monk_ind);
load(fullfile(project_dir, output_database, fname_), 'temp_input');

fname_ = sprintf('spectrogram_change_t_%d_%d_%d_%d_monkey_%d).mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2), monk_ind);
load(fullfile(project_dir, output_database, fname_), 'ttest_temp');

%%
%   Visualize F-satatistics
image(reshape([tbl_temp{4,6,:}],size(temp_input{1}, 2), size(temp_input{1}, 3)), 'CDataMapping', 'scaled')
%   Connectivity labeling
[~, temp_n]= bwlabel(reshape([tbl_temp{3,7,:}],size(temp_input{1}, 2), size(temp_input{1}, 3)) < 0.05);
%%  Mean value plot
figure()
hold on
plot_f_range = 5:9;
% plot_f_range = 8:18
plot(squeeze(mean(mean(temp_input{1}(:, plot_f_range, :), 1), 2)));
plot(squeeze(mean(mean(temp_input{2}(:, plot_f_range, :), 1), 2)));
plot(squeeze(mean(mean(temp_input{3}(:, plot_f_range, :), 1), 2)));
plot(squeeze(mean(mean(temp_input{4}(:, plot_f_range, :), 1), 2)));
legend()
%%
plot_x = [-1, 2];
plot_y = [2, 100];
%%  AS: Adult vs. Young
% temp_crit4 = reshape([tbl_temp{4,7,:}],size(temp_input{1}, 2), size(temp_input{1}, 3))<0.05;
temp_crit2 = squeeze((mean(temp_input{4}, 1) - mean(temp_input{2}, 1)) > 0);
temp_crit3 = squeeze((mean(temp_input{4}, 1) - mean(temp_input{2}, 1)) < 0);
temp_crit1 = reshape(ttest_temp(:, 2),size(temp_input{1}, 2), size(temp_input{1}, 3))<0.05;
temp_value = squeeze((mean(temp_input{4}, 1) - mean(temp_input{2}, 1)));
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_crit2 - temp_crit1.*temp_crit3, [-1, 1])
zw_set_plot_5_5(1);
title('AS')
fig_name = sprintf('as_val_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')
% 
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_value, [-.2, .2]);
zw_set_plot_5_5(4);
title('AS')
fig_name = sprintf('as_t_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')

% figure()
% image(temp_crit1, 'CDataMapping', 'scaled')
% figure()
% image(temp_crit2, 'CDataMapping', 'scaled')
% figure()
%%  PS: Adult vs. Young
% temp_crit4 = reshape([tbl_temp{4,7,:}],size(temp_input{1}, 2), size(temp_input{1}, 3))<0.05;
temp_crit2 = squeeze((mean(temp_input{3}, 1) - mean(temp_input{1}, 1)) > 0);
temp_crit3 = squeeze((mean(temp_input{3}, 1) - mean(temp_input{1}, 1)) < 0);
temp_crit1 = reshape(ttest_temp(:, 1),size(temp_input{1}, 2), size(temp_input{1}, 3))<0.05;
temp_value = squeeze((mean(temp_input{3}, 1) - mean(temp_input{1}, 1)));
% 
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_crit2 - temp_crit1.*temp_crit3, [-1, 1])
zw_set_plot_5_5(1);
title('PS')
fig_name = sprintf('ps_val_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')
% 
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_value, [-.2, .2]);
zw_set_plot_5_5(4);
title('PS')
fig_name = sprintf('ps_t_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')
% figure()
% image(temp_crit1, 'CDataMapping', 'scaled')
% figure()
% image(temp_crit2, 'CDataMapping', 'scaled')
% figure()
%%  Main effect of Stage
fig_size = [5 1 8 5];
tbl_order = 2;
temp_value = reshape([tbl_temp{tbl_order,6,:}],size(temp_input{1}, 2), size(temp_input{1}, 3));
temp_crit2 = squeeze((mean(cat(1, temp_input{2}, temp_input{4}), 1) - mean(cat(1, temp_input{1}, temp_input{3}), 1))>0);
temp_crit3 = squeeze((mean(cat(1, temp_input{2}, temp_input{4}), 1) - mean(cat(1, temp_input{1}, temp_input{3}), 1))<0);
temp_crit1 = reshape([tbl_temp{tbl_order,7,:}],size(temp_input{1}, 2), size(temp_input{1}, 3))<0.05;
% 
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_crit2 - temp_crit1.*temp_crit3, [-1,1])
zw_set_plot_5_5(1)
fig_name = sprintf('stage_f_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')
%
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_crit2.*temp_value, [0 25])
zw_set_plot_5_5(3)
fig_name = sprintf('stage_f_p_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')
% 
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_crit3.*temp_value, [0 25])
zw_set_plot_5_5(3)
fig_name = sprintf('stage_f_n_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')
%%  Main effect of Task
fig_size = [5 1 8 5];
tbl_order = 3;
temp_value = reshape([tbl_temp{tbl_order,6,:}],size(temp_input{1}, 2), size(temp_input{1}, 3));
temp_crit2 = squeeze((mean(cat(1, temp_input{2}, temp_input{4}), 1) - mean(cat(1, temp_input{1}, temp_input{3}), 1))>0);
temp_crit3 = squeeze((mean(cat(1, temp_input{2}, temp_input{4}), 1) - mean(cat(1, temp_input{1}, temp_input{3}), 1))<0);
temp_crit1 = reshape([tbl_temp{tbl_order,7,:}],size(temp_input{1}, 2), size(temp_input{1}, 3))<0.05;
% 
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_crit2 - temp_crit1.*temp_crit3, [-1,1])
zw_set_plot_5_5(2)
fig_name = sprintf('task_f_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')
%
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_crit2.*temp_value, [0 25])
zw_set_plot_5_5(3)
fig_name = sprintf('task_f_p_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')
% 
figure('Units', 'inches','Position', fig_size);
imagesc(plot_x, plot_y, temp_crit1.*temp_crit3.*temp_value, [0 25])
zw_set_plot_5_5(3)
fig_name = sprintf('task_f_n_monkey_%d', monk_ind);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300')