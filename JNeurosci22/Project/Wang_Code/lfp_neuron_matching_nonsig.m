fname_ = 'lfp_tbl.mat';
load(fullfile(project_dir, output_database, fname_), 'lfp_tbl');
%   Load manual rejections
fname_ = 'maunal_session_rejection_flag.mat';
load(fullfile(project_dir, output_database, fname_), 'session_flag_list', 'target_cwts');
%%  Converts maunally rejected session flag to 0 and 1
valid_lfp = zeros(size(lfp_tbl, 1), 1);
for i = 1:size(lfp_tbl, 1)
    if any(target_cwts(~logical(session_flag_list)) == i)
        valid_lfp(i) = 1;
    end
end     
%%  Base 2D matrix mapping LFP and neuron by electrodes, eliminating invalid sessions
mapping_mat_nonsig = zeros(size(lfp_tbl, 1), size(neuron_tbl, 1));
for i = 1:size(lfp_tbl, 1)
    i
    mapping_mat_nonsig(i, :) = prod(table2array(lfp_tbl(i, 2:5)) == table2array(neuron_tbl(:, [6, 7, 8, 3])), 2).*neuron_tbl.PFC.*~responsive;
end
mapping_mat_nonsig = mapping_mat_nonsig.*valid_lfp.*(lfp_tbl.ps_file == 1);
%%  Indexing each stage-by-animal in 't' (LFP) and 'n' (neuron)
t_nonsig = cell(2, 4);
n_nonsig = cell(2, 4);
for i_stage = 1:2
    for i_monk = 1:4
        t_nonsig{i_stage, i_monk} = find(any((lfp_tbl.stage == i_stage).*(lfp_tbl.monkey_id == i_monk).*mapping_mat_nonsig, 2))';
        n_nonsig{i_stage, i_monk} = find(any((lfp_tbl.stage == i_stage).*(lfp_tbl.monkey_id == i_monk).*mapping_mat_nonsig, 1));
    end
end
%%  Neuron tuning categorization (gamma X spiking).
temp_counts_ = {};
neuron_tuning_cat_nonsig = cell(2, 4, 2, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, gamma tuning, spiking tuning
        neuron_tuning_cat_nonsig{i_stage, i_monk, 1, 1} = intersect(n_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*~tfr_delay_tuning_flag(:, 3).*~neuron_epoch_delay_tuning_flag, 1)));        
        neuron_tuning_cat_nonsig{i_stage, i_monk, 1, 2} = intersect(n_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*~tfr_delay_tuning_flag(:, 3).*neuron_epoch_delay_tuning_flag, 1)));
        neuron_tuning_cat_nonsig{i_stage, i_monk, 2, 1} = intersect(n_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*tfr_delay_tuning_flag(:, 3).*~neuron_epoch_delay_tuning_flag, 1)));
        neuron_tuning_cat_nonsig{i_stage, i_monk, 2, 2} = intersect(n_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*tfr_delay_tuning_flag(:, 3).*neuron_epoch_delay_tuning_flag, 1)));
        temp_counts_{i_stage, i_monk} = [numel(neuron_tuning_cat_nonsig{i_stage, i_monk, 1, 1}), numel(neuron_tuning_cat_nonsig{i_stage, i_monk, 1, 2}); ...
            numel(neuron_tuning_cat_nonsig{i_stage, i_monk, 2, 1}), numel(neuron_tuning_cat_nonsig{i_stage, i_monk, 2, 2})];
    end
end
%%  Site tuning categorization (gamma X spiking)
site_tuning_cat_nonsig = cell(2, 4, 2, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, gamma tuning, spiking tuning
        site_tuning_cat_nonsig{i_stage, i_monk, 1, 1} = intersect(t_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*~tfr_delay_tuning_flag(:, 3).*~neuron_epoch_delay_tuning_flag, 2))');        
        site_tuning_cat_nonsig{i_stage, i_monk, 1, 2} = intersect(t_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*~tfr_delay_tuning_flag(:, 3).*neuron_epoch_delay_tuning_flag, 2))');
        site_tuning_cat_nonsig{i_stage, i_monk, 2, 1} = intersect(t_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*tfr_delay_tuning_flag(:, 3).*~neuron_epoch_delay_tuning_flag, 2))');
        site_tuning_cat_nonsig{i_stage, i_monk, 2, 2} = intersect(t_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*tfr_delay_tuning_flag(:, 3).*neuron_epoch_delay_tuning_flag, 2))');
    end
end
 
%%  Neuron Gamma power modulation categorization.
neuron_mod_cat_nonsig = cell(2, 4, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey,   delay gamma mod
        neuron_mod_cat_nonsig{i_stage, i_monk, 1} = intersect(n_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*~tfr_delay_modulation_flag(:, 3, 1), 1)));
        neuron_mod_cat_nonsig{i_stage, i_monk, 2} = intersect(n_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*tfr_delay_modulation_flag(:, 3, 1), 1)));
    end
end
%%  Site Gamma power modulation categorization.
site_mod_cat_nonsig = cell(2, 4, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, delay gamma mod
        site_mod_cat_nonsig{i_stage, i_monk, 1} = intersect(t_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*~tfr_delay_modulation_flag(:, 3, 1), 2))');
        site_mod_cat_nonsig{i_stage, i_monk, 2} = intersect(t_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*tfr_delay_modulation_flag(:, 3, 1), 2))');
        site_mod_cat_nonsig{i_stage, i_monk, 3} = intersect(t_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*~tfr_delay_modulation_flag(:, 4, 1), 2))');
        site_mod_cat_nonsig{i_stage, i_monk, 4} = intersect(t_nonsig{i_stage, i_monk}, find(any(mapping_mat_nonsig.*tfr_delay_modulation_flag(:, 4, 1), 2))');        
    end
end
%%
save(fullfile(project_dir, output_database, 'neuron_site_categories_nonsig.mat'), 'site_mod_cat_nonsig', 'site_tuning_cat_nonsig', 'neuron_mod_cat_nonsig', 'neuron_tuning_cat_nonsig', 'mapping_mat_nonsig');
