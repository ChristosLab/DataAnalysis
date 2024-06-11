clear
zw_setpath
%%
fname_ = 'cwt_error_repo.mat';
fname_ = 'cwt_error_repo_edited.mat';
load(fullfile(project_dir, output_database, fname_), 'cwt_error_repo');
fname_ = 'complete_cwt_repo_3_2_2-21.mat';
load(fullfile(project_dir, output_database, fname_), 'cwt_repo');
fname_ = 'lfp_tbl.mat';
load(fullfile(project_dir, output_database, fname_), 'lfp_tbl');
fname_ = 'neuron_site_categories.mat';
load(fullfile(project_dir, output_database, fname_), 'site_mod_cat', 'site_tuning_cat', 'neuron_mod_cat', 'neuron_tuning_cat', 'mapping_mat', 'responsive');
fname_ = 'lfp_neuron_matching.mat';
load(fullfile(project_dir, output_database, fname_));
fname_ = 'maunal_session_rejection_flag_error.mat';
load(fullfile(project_dir, output_database, fname_), 'session_flag_list');
%%  CWT parameters
fs = 500;
cue_dur = [-1 + 1/fs, 4]; % Signal length 
sac_dur = [-2 + 1/fs, 3];
normalizer = 4;
kernel_flag = 3;
down_sample = 10;
f_range = 2:2:128;
avg_method = 1; %   0 - arithmetic; 1 - geometric
cue_n_sample = (diff(cue_dur)*fs + 1)/down_sample;
sac_n_sample = (diff(sac_dur)*fs + 1)/down_sample;
target_frs = [4,8; 8,16; 16, 32; 32, numel(f_range)];
baseline_bin = [26, 50];
%%
t_y = [t{1, :}];
t_a = [t{2, :}];
%%
t_ya = {t_y, t_a};
error_matched_cwt = cell(size(t_ya));
error_seesion_num = [];
for i = 1:numel(t_ya)
    session_counter = 0;
    temp_cwt_match  = zeros(0, numel(f_range), cue_n_sample);
    temp_cwt_error  = zeros(0, numel(f_range), cue_n_sample);
    for j = 1:numel(t_ya{i})
        current_session = t_ya{i}(j);
        if isempty(cwt_error_repo(current_session).class(1).norm) && numel(cwt_error_repo(current_session).class) == 1
            continue
        end
        if session_flag_list(current_session)
            continue
        end
        error_seesion_num = [error_seesion_num, current_session];
        session_counter = session_counter + 1;
        single_cwt_cue        = zeros(0, numel(f_range), cue_n_sample);
        single_cwt_cue_error  = zeros(0, numel(f_range), cue_n_sample);
        for k = 1:numel(cwt_error_repo(current_session).class)
            if ~isempty(cwt_error_repo(current_session).class(k).cue_cwt)
                single_cwt_cue        = [single_cwt_cue; cwt_repo(current_session).class(k).cue_cwt./mean(cwt_repo(current_session).class(k).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
                single_cwt_cue_error  = [single_cwt_cue_error; cwt_error_repo(current_session).class(k).cue_cwt./mean(cwt_error_repo(current_session).class(k).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
            end
        end
        temp_cwt_match = [temp_cwt_match; nanmean(single_cwt_cue, 1)];
        temp_cwt_error = [temp_cwt_error; nanmean(single_cwt_cue_error, 1)];
    end
    error_matched_cwt{i} = {temp_cwt_match, temp_cwt_error};
end
%%
for i = 1:size(target_frs, 1)
    lfp_mod_cue_match_y(:, :, i) = squeeze(nanmean(error_matched_cwt{1}{1}(:, target_frs(i, 1):target_frs(i, 2), :), 2));
    lfp_mod_cue_error_y(:, :, i) = squeeze(nanmean(error_matched_cwt{1}{2}(:, target_frs(i, 1):target_frs(i, 2), :), 2));
    lfp_mod_cue_match_a(:, :, i) = squeeze(nanmean(error_matched_cwt{2}{1}(:, target_frs(i, 1):target_frs(i, 2), :), 2));
    lfp_mod_cue_error_a(:, :, i) = squeeze(nanmean(error_matched_cwt{2}{2}(:, target_frs(i, 1):target_frs(i, 2), :), 2));
end
%%
fname_ = 'lfp_mod_cue_error';
save(fullfile(project_dir, output_database, fname_), 'lfp_mod_cue_match_y', 'lfp_mod_cue_error_y', 'lfp_mod_cue_match_a', 'lfp_mod_cue_error_a');
%%
fname_ = 'lfp_mod_cue_error';
load(fullfile(project_dir, output_database, fname_), 'lfp_mod_cue_match_y', 'lfp_mod_cue_error_y', 'lfp_mod_cue_match_a', 'lfp_mod_cue_error_a');
%%
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
for band_i_ = 1:4
    err_mod_plot({lfp_mod_cue_match_y(:, :, band_i_), ...
        lfp_mod_cue_error_y(:, :, band_i_), ...
        lfp_mod_cue_match_a(:, :, band_i_), ...
        lfp_mod_cue_error_a(:, :, band_i_)}, ...
        linspace(-1,4,250), ...
        {'Ado c', 'Ado e', 'Adu c', 'Adu e'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, [0.7, 2.1])
    %     {'Adolescent preferred', 'Adolescent diametric', 'Adult preferred', 'Adult diametric'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
    title(title_st{band_i_})
    xlabel('Time from cue onset (s)')
    % ylim([0.7, 2.4])
    xlim([-1, 2.5])
    if band_i_ ~= 2
        legend off
    end
    fig_name = ['lfp_pow_band_', num2str(band_i_), 'power_modulation_error_match', '_draft'];
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
    savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
end
%%