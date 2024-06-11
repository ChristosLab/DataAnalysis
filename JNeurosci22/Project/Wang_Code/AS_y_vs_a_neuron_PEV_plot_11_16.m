%%  Young vs. adult AS cue-aligned
inputs = {neuron_pev_cue(find(needed_neuron.*(neuron_tbl.task_id == 2).*(neuron_tbl.aged == 0)), :), ...
    neuron_pev_cue(find(needed_neuron.*(neuron_tbl.task_id == 2).*(neuron_tbl.aged == 1)), :)};
pev_plot(inputs, t_bins_cue, {'Adolescent', 'Adult'}, {'r', 'b'})
xlabel('Time from cue onset (s)')
title('PEV of neurons at different stages (anti-saccade)')
fig_name = 'pev_stage_as_cue_aligned';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Young vs. adult AS sac-aligned
inputs = {neuron_pev_sac(find(needed_neuron.*(neuron_tbl.task_id == 2).*(neuron_tbl.aged == 0)), :), ...
    neuron_pev_sac(find(needed_neuron.*(neuron_tbl.task_id == 2).*(neuron_tbl.aged == 1)), :)};
pev_plot(inputs, t_bins_sac, {'Adolescent', 'Adult'}, {'r', 'b'})
xlabel('Time from saccadic onset (s)')
title('PEV of neurons at different stages (anti-saccade)')
fig_name = 'pev_stage_as_sac_aligned';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');