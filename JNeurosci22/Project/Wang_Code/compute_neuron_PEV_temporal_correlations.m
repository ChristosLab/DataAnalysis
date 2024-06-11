%%  Compute Pearson correlation between neuron PEV and LFP modulation of each band
%                             neuron                           X band 
neuron_pev_lfp_mod_cor_r = nan(size(neuron_pev_cue, 1), size(lfp_mod_cue, 3));
neuron_pev_lfp_mod_cor_p = nan(size(neuron_pev_cue, 1), size(lfp_mod_cue, 3));

cor_window_ = 76:150;
for i = 1:size(mapping_mat, 2) %    looping through
    if any(mapping_mat(:, i))
        for j = 1:size(lfp_mod_cue, 3)
%             [cor_r_, cor_p_] = corrcoef(lfp_mod_cue(find(mapping_mat(:, i)), cor_window_, j), neuron_pev_cue(i, cor_window_));
%             neuron_pev_lfp_mod_cor_r(i, j) = cor_r_(2, 1);
%             neuron_pev_lfp_mod_cor_p(i, j) = cor_p_(2, 1);
            [neuron_pev_lfp_mod_cor_r(i, j), neuron_pev_lfp_mod_cor_p(i, j)] = corr(lfp_mod_cue(find(mapping_mat(:, i)), cor_window_, j)', neuron_pev_cue(i, cor_window_)', 'type', 'Spearman');
            if isnan(neuron_pev_lfp_mod_cor_r(i, j))
                neuron_pev_lfp_mod_cor_r(i, j) = 0;
            end
        end
    end
end
%%  Compute Pearson correlation between neuron PEV and LFP PEV of each band
%                             neuron                           X band 
neuron_pev_lfp_pev_cor_r = nan(size(neuron_pev_cue, 1), size(lfp_pev_cue, 3));
neuron_pev_lfp_pev_cor_p = nan(size(neuron_pev_cue, 1), size(lfp_pev_cue, 3));

cor_window_ = 76:150;
for i = 1:size(mapping_mat, 2) %    looping through
    if any(mapping_mat(:, i))
        for j = 1:size(lfp_mod_cue, 3)
            [cor_r_, cor_p_] = corrcoef(lfp_pev_cue(find(mapping_mat(:, i)), cor_window_, j), neuron_pev_cue(i, cor_window_));
            neuron_pev_lfp_pev_cor_r(i, j) = cor_r_(2, 1);
            neuron_pev_lfp_pev_cor_p(i, j) = cor_p_(2, 1);
        end
    end
end
%%  Resample the population to contruct a significance criterion (NOT FINISHED)
% neuron_pop_for_cor = [neuron_tuning_cat{:, [1,3,4], :, :}];
% n_resample = 50000;
% cor_window_ = 76:150;
% resampled_neuron_pev_lfp_pev_cor_r = zeros(size(lfp_mod_cue, 3), n_resample);
% resampled_neuron_pev_lfp_pev_cor_p = zeros(size(lfp_mod_cue, 3), n_resample);
% 
% for i = 1:n_resample %    looping through
%         for j = 1:size(lfp_mod_cue, 3)
%             [cor_r_, cor_p_] = corrcoef(lfp_pev_cue(, cor_window_, j), neuron_pev_cue(, cor_window_));
%             resampled_neuron_pev_lfp_pev_cor_r(i, j) = cor_r_(2, 1);
%             resampled_neuron_pev_lfp_pev_cor_p(i, j) = cor_p_(2, 1);
%         end
%     end
% end
%%
numel(intersect([neuron_tuning_cat{1, [1,3,4], :, :}], find((neuron_pev_lfp_mod_cor_p(:, 3) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 3) < 0)), 1))

numel(intersect([neuron_tuning_cat{1, [1,3,4], :, :}], find((neuron_pev_lfp_mod_cor_p(:, 3) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 3) < 0)), 2))
numel(intersect([neuron_tuning_cat{1, [1,3,4], :, :}], find((neuron_pev_lfp_mod_cor_p(:, 3) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 3) < 0)), 3))
numel(intersect([neuron_tuning_cat{1, [1,3,4], :, :}], find((neuron_pev_lfp_mod_cor_p(:, 3) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 3) < 0)), 4))
%%
figure; histogram(neuron_pev_lfp_mod_cor_r(intersect([neuron_tuning_cat{1, [1,3,4], :, 2}], find(neuron_pev_lfp_mod_cor_p(:, 3) < 0.05)), 3))
%%
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 3), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, 1}], 3), -1:0.1:1, 'Normalization', 'probability')
%%
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 3), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, 1}], 3), -1:0.1:1, 'Normalization', 'probability')
%%
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 3), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 3), -1:0.1:1, 'Normalization', 'probability')
%%  PEV_VS_MOD: All adolescent vs. all adult
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, :}], 1), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, :}], 1), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, :}], 2), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, :}], 2), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, :}], 3), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, :}], 3), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, :}], 4), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, :}], 4), -1:0.1:1, 'Normalization', 'probability')
%%  PEV_VS_MOD: Informative adolescent vs. informative adult
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 1), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 1), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 2), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 2), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 3), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 3), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 4), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 4), -1:0.1:1, 'Normalization', 'probability')
%%  PEV_VS_PEV: All adolescent vs. all adult
figure; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{1, [1,3,4], :, :}], 1), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, :}], 1), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{1, [1,3,4], :, :}], 2), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, :}], 2), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{1, [1,3,4], :, :}], 3), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, :}], 3), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{1, [1,3,4], :, :}], 4), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, :}], 4), -1:0.1:1, 'Normalization', 'probability')
%%  PEV_VS_PEV: Informative adolescent vs. informative adult
figure; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 1), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 1), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 2), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 2), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 3), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 3), -1:0.1:1, 'Normalization', 'probability')
figure; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 4), -1:0.1:1, 'Normalization', 'probability')
hold on; histogram(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 4), -1:0.1:1, 'Normalization', 'probability')
%%  Scatter: neuron PEV and corr:pev_pev
figure
plot(mean(neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 2}], cor_window_), 2), neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 4), '*')
%%
figure
plot(mean(neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 2}], cor_window_), 2), neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 4), '*')
%%
figure
plot(mean(best_psth_raw([neuron_tuning_cat{2, [1,3,4], :, 2}], 31:60), 2), neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 4), '*')
%%
figure
plot(max(neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 2}], 31:60), 2), mean(mean(temp_cwt([neuron_tuning_cat{2, [1,3,4], :, 2}], target_frs(3,1):target_frs(3,2), 76:150), 3), 2), '*')
%%  Max neuron PEV to Mean Gamma power corr (Adolescent)
figure;
plot(max(neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, :}], 31:60), [], 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{1, [1,3,4], :, :}]), 76:150, 3), 2), '*')
%%  Max neuron PEV to Mean Gamma power corr (Adult)
figure;
plot(max(neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, :}], 31:60), [], 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{2, [1,3,4], :, :}]), 76:150, 3), 2), '*')
%%  Average neuron PEV to Mean Gamma power corr (Adolescent)
figure;
plot(mean(neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, :}], 31:60), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{1, [1,3,4], :, :}]), 76:150, 3), 2), '*')
%%  Average neuron PEV to Mean Gamma power corr (Adult)
figure;
plot(mean(neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, :}], 31:60), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{2, [1,3,4], :, :}]), 76:150, 3), 2), '*')
%%
hist_plot({neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 3), neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 3)}, ...
    -1:0.1:1, {'Adolescent', 'Adult'}, {'b', 'r'}, 'probability');
title('Neuron PEV ~ gamma power')
fig_name = 'neuron_pev_gamma_mod_hist_draft'
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%
hist_plot({neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, [1,3,4], :, 2}], 4), neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, [1,3,4], :, 2}], 4)}, ...
    -1:0.1:1, {'Adolescent', 'Adult'}, {'b', 'r'}, 'probability');
title('Neuron PEV ~ high-gamma power')
fig_name = 'neuron_pev_high_gamma_mod_hist_draft'
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%
numel(intersect([neuron_tuning_cat{1, :, :, 2}], find((neuron_pev_lfp_mod_cor_p(:, 3) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 3) < 0))))
numel(intersect([neuron_tuning_cat{1, :, :, 2}], find((neuron_pev_lfp_mod_cor_p(:, 3) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 3) > 0))))
numel(intersect([neuron_tuning_cat{2, :, :, 2}], find((neuron_pev_lfp_mod_cor_p(:, 3) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 3) < 0))))
numel(intersect([neuron_tuning_cat{2, :, :, 2}], find((neuron_pev_lfp_mod_cor_p(:, 3) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 3) > 0))))
%%
numel(intersect([neuron_tuning_cat{1, :, :, 2}], find((neuron_pev_lfp_mod_cor_p(:, 4) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 4) < 0))))
numel(intersect([neuron_tuning_cat{1, :, :, 2}], find((neuron_pev_lfp_mod_cor_p(:, 4) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 4) > 0))))
numel(intersect([neuron_tuning_cat{2, :, :, 2}], find((neuron_pev_lfp_mod_cor_p(:, 4) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 4) < 0))))
numel(intersect([neuron_tuning_cat{2, :, :, 2}], find((neuron_pev_lfp_mod_cor_p(:, 4) < 0.05).*(neuron_pev_lfp_mod_cor_r(:, 4) > 0))))
%%
clc
[h, p, ~, stats] = ttest2(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, :, :, 2}], 3), neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, :, :, 2}], 3))
%%
clc
[h, p, ~, stats] = ttest2(neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{1, :, :, 2}], 4), neuron_pev_lfp_mod_cor_r([neuron_tuning_cat{2, :, :, 2}], 4))