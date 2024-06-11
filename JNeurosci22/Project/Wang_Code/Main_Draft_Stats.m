cue_dur = 51:75;
delay_dur = 76:150;
cue_dur = delay_dur;
%%  Band power mean and stage difference
disp('Apha:')
mean(mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], cue_dur, 1))) - 1
mean(mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], cue_dur, 1))) - 1
disp('Beta:')
mean(mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], cue_dur, 2))) - 1
mean(mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], cue_dur, 2))) - 1
disp('Gamma:')
mean(mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], cue_dur, 3))) - 1
mean(mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], cue_dur, 3))) - 1
disp('High-gamma:')
mean(mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], cue_dur, 4))) - 1
mean(mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], cue_dur, 4))) - 1
%%  Test each stage effect
cue_dur = 51:75;
delay_dur = 76:150;
% cue_dur = delay_dur;
% 
disp('Apha:')
ts_1 = mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], cue_dur, 1), 2) - 1;
ts_2 = mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], cue_dur, 1), 2) - 1;
[h_, p_, ~, stats] = ttest2(ts_1, ts_2)
disp('Beta:')
ts_1 = mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], cue_dur, 2), 2) - 1;
ts_2 = mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], cue_dur, 2), 2) - 1;
[h_, p_, ~, stats] = ttest2(ts_1, ts_2)
disp('Gamma:')
ts_1 = mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], cue_dur, 3), 2) - 1;
ts_2 = mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], cue_dur, 3), 2) - 1;
[h_, p_, ~, stats] = ttest2(ts_1, ts_2)
disp('High-gamma:')
ts_1 = mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], cue_dur, 4), 2) - 1;
ts_2 = mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], cue_dur, 4), 2) - 1;
[h_, p_, ~, stats] = ttest2(ts_1, ts_2)
%%
clc
disp('Apha:')
ts_1 = mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], delay_dur, 1), 2) - 1;
ts_2 = mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], delay_dur, 1), 2) - 1;
[h_, p_, ~, stats] = ttest2(ts_1, ts_2)
disp('Beta:')
ts_1 = mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], delay_dur, 2), 2) - 1;
ts_2 = mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], delay_dur, 2), 2) - 1;
[h_, p_, ~, stats] = ttest2(ts_1, ts_2)
disp('Gamma:')
ts_1 = mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], delay_dur, 3), 2) - 1;
ts_2 = mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], delay_dur, 3), 2) - 1;
[h_, p_, ~, stats] = ttest2(ts_1, ts_2)
disp('High-gamma:')
ts_1 = mean(lfp_mod_cue([site_mod_cat{1, [1,3 , 4], 1:2}], delay_dur, 4), 2) - 1;
ts_2 = mean(lfp_mod_cue([site_mod_cat{2, [1,3 , 4], 1:2}], delay_dur, 4), 2) - 1;
[h_, p_, ~, stats] = ttest2(ts_1, ts_2)
%%  Gamma power comparison at gamma modulated sites between adolescent and adult
set_gamma_power_mod_adolescent = mean(temp_cwt_b([site_mod_cat{1, [1,3 , 4], 2}], 76:150), 2);
set_gamma_power_mod_adult = mean(temp_cwt_b([site_mod_cat{2, [1,3 , 4], 2}], 76:150), 2);
[h_, p_, ~, stats] = ttest2(set_gamma_power_mod_adolescent, set_gamma_power_mod_adult)
% 
%%  Modulated site counting
disp('young gamma mod')
numel([site_mod_cat{1, [1,3,4], 2}])
disp('young gamma total')
numel([site_mod_cat{1, [1,3,4], 1:2}])
disp('adult gamma mod')
numel([site_mod_cat{2, [1,3,4], 2}])
disp('adult gamma total')
numel([site_mod_cat{2, [1,3,4], 1:2}])
disp('young h gamma mod')
numel([site_mod_cat{1, [1,3,4], 4}])
disp('young h gamma total')
numel([site_mod_cat{1, [1,3,4], 3:4}])
disp('adult h gamma mod')
numel([site_mod_cat{2, [1,3,4], 4}])
disp('adult h gamma total')
numel([site_mod_cat{2, [1,3,4], 3:4}])
%%  Modulated delay power compare
clc
band_i = 3;
set_gamma_power_mod_adolescent = mean(lfp_mod_cue([site_mod_cat{1, [1, 3, 4], 2}], 76:150, band_i), 2);
set_gamma_power_mod_adult = mean(lfp_mod_cue([site_mod_cat{2, [1, 3, 4], 2}], 76:150, band_i), 2);
[h_, p_, ~, stats] = ttest2(set_gamma_power_mod_adolescent, set_gamma_power_mod_adult)
% 
band_i = 4;
set_gamma_power_mod_adolescent = mean(lfp_mod_cue([site_mod_cat{1, [1, 3, 4], 4}], 76:150, band_i), 2);
set_gamma_power_mod_adult = mean(lfp_mod_cue([site_mod_cat{2, [1, 3, 4], 4}], 76:150, band_i), 2);
[h_, p_, ~, stats] = ttest2(set_gamma_power_mod_adolescent, set_gamma_power_mod_adult)
%%  Spatially-tune site counting
disp('young gamma tune')
numel([site_tuning_cat{1, [1,3,4], 2, :}])
disp('young gamma total')
numel([site_tuning_cat{1, [1,3,4], 1:2, :}])
disp('adult gamma tune')
numel([site_tuning_cat{2, [1,3,4], 2, :}])
disp('adult gamma total')
numel([site_tuning_cat{2, [1,3,4], 1:2, :}])
disp('young high gamma tune')
numel([site_tuning_cat_hg{1, [1,3,4], 2, :}])
disp('young high gamma total')
numel([site_tuning_cat_hg{1, [1,3,4], 1:2, :}])
disp('adult high gamma tune')
numel([site_tuning_cat_hg{2, [1,3,4], 2, :}])
disp('adult high gamma total')
numel([site_tuning_cat_hg{2, [1,3,4], 1:2, :}])
%%  Spatially-tuned delay power compare
clc
band_i = 3;
set_gamma_power_mod_adolescent = mean(lfp_mod_cue([site_tuning_cat{1, [1, 3, 4], 2, :}], 76:150, band_i), 2);
set_gamma_power_mod_adult = mean(lfp_mod_cue([site_tuning_cat{2, [1, 3, 4], 2, :}], 76:150, band_i), 2);
[h_, p_, ~, stats] = ttest2(set_gamma_power_mod_adolescent, set_gamma_power_mod_adult)
%%
band_i = 4;
set_gamma_power_mod_adolescent = mean(lfp_mod_cue([site_tuning_cat_hg{1, [1, 3, 4], 2, :}], 76:150, band_i), 2);
set_gamma_power_mod_adult = mean(lfp_mod_cue([site_tuning_cat_hg{2, [1, 3, 4], 2, :}], 76:150, band_i), 2);
[h_, p_, ~, stats] = ttest2(set_gamma_power_mod_adolescent, set_gamma_power_mod_adult)
%%  Non spatially-tuned delay power compare
clc
band_i = 3;
set_gamma_power_mod_adolescent = mean(lfp_mod_cue([site_tuning_cat{1, [1, 3, 4], 1, :}], 76:150, band_i), 2);
set_gamma_power_mod_adult = mean(lfp_mod_cue([site_tuning_cat{2, [1, 3, 4], 1, :}], 76:150, band_i), 2);
[h_, p_, ~, stats] = ttest2(set_gamma_power_mod_adolescent, set_gamma_power_mod_adult)
% 
band_i = 4;
set_gamma_power_mod_adolescent = mean(lfp_mod_cue([site_tuning_cat_hg{1, [1, 3, 4], 1, :}], 76:150, band_i), 2);
set_gamma_power_mod_adult = mean(lfp_mod_cue([site_tuning_cat_hg{2, [1, 3, 4], 1, :}], 76:150, band_i), 2);
[h_, p_, ~, stats] = ttest2(set_gamma_power_mod_adolescent, set_gamma_power_mod_adult)

%%  Number of adolescent neuron with significant neuron_pev vs. lfp_mod correlation
i_band = 3;
i_stage = 1;
sum((neuron_pev_lfp_pev_cor_p([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band) < 0.05).*...
    (neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band)) < 0)
sum((neuron_pev_lfp_pev_cor_p([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band) < 0.05).*...
    (neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band)) > 0)
i_band = 4;
sum((neuron_pev_lfp_pev_cor_p([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band) < 0.05).*...
    (neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band)) < 0)
sum((neuron_pev_lfp_pev_cor_p([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band) < 0.05).*...
    (neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band)) > 0)
%%  Number of adult neuron with significant neuron_pev vs. lfp_mod correlation
i_band = 3;
i_stage = 2;
sum((neuron_pev_lfp_pev_cor_p([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band) < 0.05).*...
    (neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band)) < 0)
sum((neuron_pev_lfp_pev_cor_p([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band) < 0.05).*...
    (neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band)) > 0)
i_band = 4;
sum((neuron_pev_lfp_pev_cor_p([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band) < 0.05).*...
    (neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band)) < 0)
sum((neuron_pev_lfp_pev_cor_p([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band) < 0.05).*...
    (neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band)) > 0)
%%  Number of informative neurons with gamma mod
stage_i = 1;
numel(intersect([neuron_mod_cat{stage_i, [1,3,4], 2}], [neuron_tuning_cat{stage_i,[1,3,4], :, 2}]))
numel([neuron_tuning_cat{stage_i,[1,3,4], :, 2}])
numel(intersect([neuron_mod_cat{stage_i, [1,3,4], 2}], [neuron_tuning_cat{stage_i,[1,3,4], :, 1}]))
numel([neuron_tuning_cat{stage_i,[1,3,4], :, 1}])
%%
stage_i = 2;
numel(intersect([neuron_mod_cat{stage_i, [1,3,4], 2}], [neuron_tuning_cat{stage_i,[1,3,4], :, 2}]))
numel([neuron_tuning_cat{stage_i,[1,3,4], :, 2}])
numel(intersect([neuron_mod_cat{stage_i, [1,3,4], 2}], [neuron_tuning_cat{stage_i,[1,3,4], :, 1}]))
numel([neuron_tuning_cat{stage_i,[1,3,4], :, 1}])

%%  neuron_pev vs. lfp_mod correlation means
i_band = 3;
i_stage = 1;
mean(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band))
i_stage = 2;
mean(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band))
i_band = 4;
i_stage = 1;
mean(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band))
i_stage = 2;
mean(neuron_pev_lfp_pev_cor_r([neuron_tuning_cat{i_stage, [1,3,4], :, 2}], i_band))



%%  Permutation test for firing rate difference significance
set_psth_raw_non_mod_adolescent = mean(best_psth_raw([neuron_mod_cat{1, [1,3,4], 1}], 31:60), 2);
set_psth_raw_non_mod_adult = mean(best_psth_raw([neuron_mod_cat{2, [1,3,4], 1}], 31:60), 2);
set_psth_raw_mod_adolescent = mean(best_psth_raw([neuron_mod_cat{1, [1,3,4], 2}], 31:60), 2);
set_psth_raw_mod_adult = mean(best_psth_raw([neuron_mod_cat{2, [1,3,4], 2}], 31:60), 2);
% 
set_psth_raw_adult = mean(best_psth_raw([neuron_mod_cat{2, [1,3,4], :}], 31:60), 2);
% 
set_psth_1s_baseline_non_mod_adolescent = mean(best_psth_1s_baseline([neuron_mod_cat{1, [1,3,4], 1}], 31:60), 2);
set_psth_1s_baseline_non_mod_adult = mean(best_psth_1s_baseline([neuron_mod_cat{2, [1,3,4], 1}], 31:60), 2);
set_psth_1s_baseline_mod_adolescent = mean(best_psth_1s_baseline([neuron_mod_cat{1, [1,3,4], 2}], 31:60), 2);
set_psth_1s_baseline_mod_adult = mean(best_psth_1s_baseline([neuron_mod_cat{2, [1,3,4], 2}], 31:60), 2);
%%  Mean of firing rate differences: adolescent_gamma vs adolescent_no_gamma
f_ = zw_permute_test(set_psth_raw_mod_adolescent, set_psth_raw_non_mod_adolescent, 100000)
%%  Mean of evoked firing rate differences: adolescent_gamma vs adolescent_no_gamma
f_ = zw_permute_test(set_psth_1s_baseline_mod_adolescent, set_psth_1s_baseline_non_mod_adolescent, 100000)
%%  Mean of firing rate differences: adolescent_no_gamma vs adult
f_ = zw_permute_test(set_psth_raw_mod_adult, set_psth_raw_adult, 100000)
%%  Mean of neuronal PEV differences: adolescent_gamma vs adult_gamma
set_neuron_pev_mod_adolescent = mean(neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 2}], 31:60), 2);
set_neuron_pev_mod_adult = mean(neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 2}], 31:60), 2);
set_neuron_pev_non_mod_adolescent = mean(neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 1}], 31:60), 2);
set_neuron_pev_non_mod_adult = mean(neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 1}], 31:60), 2);
[h_, p_] = ttest2(set_neuron_pev_mod_adolescent, set_neuron_pev_mod_adult)
%%  ANOVA Neuron PEV: Stage-by-Gamma mod
anova_input_cell = {set_neuron_pev_non_mod_adolescent', set_neuron_pev_non_mod_adult'; set_neuron_pev_mod_adolescent', set_neuron_pev_mod_adult'};
zw_anova_from_cell(anova_input_cell)
%%  Mean of neuronal PEV differences: adolescent_gamma vs adult_gamma (SIG + NONSIG)
set_neuron_pev_mod_adolescent = mean(neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 2}, neuron_mod_cat_nonsig{1, [1,3,4], 2}], 31:60), 2);
set_neuron_pev_mod_adult = mean(neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 2}, neuron_mod_cat_nonsig{2, [1,3,4], 2}], 31:60), 2);
set_neuron_pev_non_mod_adolescent = mean(neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 1}, neuron_mod_cat_nonsig{1, [1,3,4], 1}], 31:60), 2);
set_neuron_pev_non_mod_adult = mean(neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 1}, neuron_mod_cat_nonsig{2, [1,3,4], 1}], 31:60), 2);
[h_, p_] = ttest2(set_neuron_pev_mod_adolescent, set_neuron_pev_mod_adult)
%%  ANOVA Neuron PEV: Stage-by-Gamma mod (SIG + NONSIG)
anova_input_cell = {set_neuron_pev_non_mod_adolescent', set_neuron_pev_non_mod_adult'; set_neuron_pev_mod_adolescent', set_neuron_pev_mod_adult'};
zw_anova_from_cell(anova_input_cell)

%%  ANOVA Firing rate (raw): Stage-by-Gamma mod
anova_input_cell = {set_psth_raw_non_mod_adolescent', set_psth_raw_non_mod_adult'; set_psth_raw_mod_adolescent', set_psth_raw_mod_adult'};
zw_anova_from_cell(anova_input_cell)
%
% 
% 
%%  ANOVA Gamma mod: Stage-by-neuron info 
band_i_ = 3
set_gamma_mod_info_adolescent = mean(lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], 76:150, band_i_), 2);
set_gamma_mod_info_adult = mean(lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], 76:150, band_i_), 2);
set_gamma_mod_non_info_adolescent = mean(lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], 76:150, band_i_), 2);
set_gamma_mod_non_info_adult = mean(lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], 76:150, band_i_), 2);
%%
anova_input_cell = {set_gamma_mod_non_info_adolescent', set_gamma_mod_non_info_adult'; set_gamma_mod_info_adolescent', set_gamma_mod_info_adult'};
zw_anova_from_cell(anova_input_cell)

%%  ANOVA Gamma PEV: Stage-by-neuron info 
band_i_ = 4
set_gamma_pev_info_adolescent = mean(lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 2}], 76:150, band_i_), 2);
set_gamma_pev_info_adult = mean(lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 2}], 76:150, band_i_), 2);
set_gamma_pev_non_info_adolescent = mean(lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 1}], 76:150, band_i_), 2);
set_gamma_pev_non_info_adult = mean(lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 1}], 76:150, band_i_), 2);
%%
anova_input_cell = {set_gamma_pev_non_info_adolescent', set_gamma_pev_non_info_adult'; set_gamma_pev_info_adolescent', set_gamma_pev_info_adult'};
zw_anova_from_cell(anova_input_cell)
%%  Spearman correlation between delay period evoked firing and high gamma power
[r_, p_] = corr(mean(best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, :}], 76:150), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{2, [1,3,4], :, :}]), 76:150, 4), 2), 'type', 'Spearman')
%%  Correlation between delay period raw firing and high gamma power
band_i = 4;
% corr_method_ = 'Pearson';
corr_method_ = 'Spearman';
[r_adolescent_ p_] = corr(mean(best_psth_raw([neuron_tuning_cat{1, [1,3,4], :, :}], 76:150), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{1, [1,3,4], :, :}]), 76:150, band_i), 2), 'type', corr_method_)
[r_adult_, p_] = corr(mean(best_psth_raw([neuron_tuning_cat{2, [1,3,4], :, :}], 76:150), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{2, [1,3,4], :, :}]), 76:150, band_i), 2), 'type', corr_method_)
disp(newline)
compare_correlation_coefficients(r_adolescent_, r_adult_, numel([neuron_tuning_cat{1, [1,3,4], :, :}]), numel([neuron_tuning_cat{2, [1,3,4], :, :}]))
disp(newline)
%%  Correlation between delay period evoked firing and gamma power
clc
band_i = 3;
% corr_method_ = 'Pearson';
corr_method_ = 'Spearman';
[r_adolescent_, p_] = corr(mean(best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, :}], 76:150), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{1, [1,3,4], :, :}]), 76:150, band_i), 2), 'type', corr_method_)
[r_adult_, p_] = corr(mean(best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, :}], 76:150), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{2, [1,3,4], :, :}]), 76:150, band_i), 2), 'type', corr_method_)
disp(newline)
compare_correlation_coefficients(r_adolescent_, r_adult_, numel([neuron_tuning_cat{1, [1,3,4], :, :}]), numel([neuron_tuning_cat{2, [1,3,4], :, :}]))
disp(newline)
%%  Correlation between delay period evoked firing and gamma power (reduced set)
clc
band_i = 4;
% corr_method_ = 'Pearson';
corr_method_ = 'Spearman';
[r_adolescent_, p_] = corr(mean(best_psth_1s_baseline(redueced_adolescent_set, 76:150), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adolescent_set), 76:150, band_i), 2), 'type', corr_method_)
[r_adult_, p_] = corr(mean(best_psth_1s_baseline(redueced_adult_set, 76:150), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adult_set), 76:150, band_i), 2), 'type', corr_method_)
disp(newline)
compare_correlation_coefficients(r_adolescent_, r_adult_, numel(redueced_adolescent_set), numel(redueced_adult_set))
disp(newline)

%%  Correlation between delay period neuron PEV and gamma power
clc
band_i = 3;
% corr_method_ = 'Pearson';
corr_method_ = 'Spearman';
[r_adolescent_, p_] = corr(mean(neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, :}], 76:150), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{1, [1,3,4], :, :}]), 76:150, band_i), 2), 'type', corr_method_)
[r_adult_, p_] = corr(mean(neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, :}], 76:150), 2),  mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{2, [1,3,4], :, :}]), 76:150, band_i), 2), 'type', corr_method_)
disp(newline)
compare_correlation_coefficients(r_adolescent_, r_adult_, numel([neuron_tuning_cat{1, [1,3,4], :, :}]), numel([neuron_tuning_cat{2, [1,3,4], :, :}]))
disp(newline)
%%  PAC comparisons
g1 = zw_pac_extract(pac_repo, [site_tuning_cat{1, [1,2,3], :, 1}], 1:10, 7:15);
g2 = zw_pac_extract(pac_repo, [site_tuning_cat{1, [1,2,3], :, 2}], 1:10, 7:15);
g3 = zw_pac_extract(pac_repo, [site_tuning_cat{2, [1,2,3], :, 1}], 1:10, 7:15);
g4 = zw_pac_extract(pac_repo, [site_tuning_cat{2, [1,2,3], :, 2}], 1:10, 7:15);

[a, sta] = zw_anova1_from_cell({g1', g2', g3', g4'})
%%
g01 = zw_pac_extract(pac_repo, [site_tuning_cat{1, [1,2,3], :, :}], 1:10, 7:15);
g02 = zw_pac_extract(pac_repo, [site_tuning_cat{2, [1,2,3], :, :}], 1:10, 7:15);
[h_, p_] = ttest2(g01, g02)
%%  Compare informative nueron PEV
[h_, p_, ~, stats_] = ttest2(squeeze(mean(neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 2}], 76:150), 2)), ...
    squeeze(mean(neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 2}], 76:150), 2)))
%%  Compare non_informative nueron PEV
[h_, p_, ~, stats_] = ttest2(squeeze(mean(neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 1}], 76:150), 2)), ...
    squeeze(mean(neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 1}], 76:150), 2)))
%%  Compare informative nueron spiking
[h_, p_, ~, stats_] = ttest2(squeeze(mean(best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 2}], 76:150), 2)), ...
    squeeze(mean(best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 2}], 76:150), 2)))