%%  Young vs. adult ODR cue-aligned (Gamma cue tuning)
fr_index = 6; % Gamma
t_lfp = (-1 + 1/epoch_fs/2):(1/epoch_fs):(2.5 - 1/epoch_fs/2) ;
single_fr_pev_ = squeeze(lfp_pev_cue(:,:,fr_index));
inputs = {...
    single_fr_pev_(find((lfp_sites_tbl.task_id == 1).*(lfp_sites_tbl.aged == 0).*(lfp_sites_tbl.responsive == 1).*(lfp_sites_tbl.n_neuron > 0)), :), ...
    single_fr_pev_(find((lfp_sites_tbl.task_id == 1).*(lfp_sites_tbl.aged == 0).*(lfp_sites_tbl.responsive == 0).*(lfp_sites_tbl.n_neuron > 0)), :), ...
    single_fr_pev_(find((lfp_sites_tbl.task_id == 1).*(lfp_sites_tbl.aged == 1).*(lfp_sites_tbl.responsive == 1).*(lfp_sites_tbl.n_neuron > 0)), :), ...
    single_fr_pev_(find((lfp_sites_tbl.task_id == 1).*(lfp_sites_tbl.aged == 1).*(lfp_sites_tbl.responsive == 0).*(lfp_sites_tbl.n_neuron > 0)), :), ...
    };
pev_plot(inputs, t_lfp, {'Adolescent w/ neuronal response', 'Adolescent w/o neuronal response', 'Adult w/ neuronal response', 'Adult w/o neuronal response'}, {'r', [0.6,0.2,0.2], 'b', [0.2, 0.2, 0.6]})
xlabel('Time from cue onset (s)')
title('PEV of LFP sites by neuronal response (ODR)')
% fig_name = 'pev_gamma_cue_odr_cue_aligned';
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%
% Notes: Sites w/ neurons show stronger, more reasonable tuning. Overall,
% there seems to be stronger delay period tuning for adolescents