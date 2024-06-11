% Read neuron epoch significance from Junda's spreadsheet
fname_ = fullfile(project_dir, output_database, 'NeuronEpoch_young.xlsx');
neuron1 = readtable(fname_, 'sheet', 'progap01');
fname_ = fullfile(project_dir, output_database, 'NeuronEpoch_adult.xlsx');
neuron2 = readtable(fname_, 'sheet', 'progap01');
fname_ = fullfile(project_dir, output_database, 'Neuron_Channels.xlsx');
neuron_channels = readtable(fname_);
% 
neuron_tbl = neuron_channels;
%
for i = 1:size(neuron_tbl, 1)
    [monkey_id_, session_id_, task_id_] = parse_fname(neuron_channels.Filename{i});
    neuron_tbl.monkey_id(i) = monkey_id_;
    neuron_tbl.session_id(i) = session_id_;
    neuron_tbl.task_id(i) = task_id_;
end
% 
new_var_range = 3:8;
for i = new_var_range
    new_var = zeros(size(neuron_tbl, 1), 1);
    neuron_tbl = addvars(neuron_tbl, new_var, 'NewVariableNames', neuron1.Properties.VariableNames(i));
end
new_var = zeros(size(neuron_tbl, 1), 1);
neuron_tbl = addvars(neuron_tbl, new_var, 'NewVariableNames', 'PFC');
for neuron_ = {neuron1, neuron2}
    for i = 1:size(neuron_{1}, 1)
        [monkey_id_, session_id_, task_id_] = parse_fname(neuron_{1}.Filename{i});
        match_ = find((neuron_tbl.Neuron == neuron_{1}.Neuron(i)).*(neuron_tbl.task_id == task_id_));
        if match_ > 0
            neuron_tbl(match_, size(neuron_channels, 2) + 4:size(neuron_channels, 2) + numel(new_var_range) + 3) = neuron_{1}(i, new_var_range);
%             neuron_tbl(match_, end - numel(new_var_range) + 1:end) = neuron_{1}(i, new_var_range);
            neuron_tbl.PFC(match_) = 1;
        end
    end
end
%%
last_young_session = [20, 26, 66, 27];
for i = 1:size(neuron_tbl, 1)
    neuron_tbl.stage(i) = 1 + (neuron_tbl.session_id(i) > last_young_session(neuron_tbl.monkey_id(i)));
end
%%
fname_ = 'neuron_tbl_w_pfc.mat';
save(fullfile(project_dir, output_database, fname_), 'neuron_tbl');
%%
function [monkey_id, session_id, task_id] = parse_fname(string_fname)
monkey_names = {'ind', 'jac', 'lem', 'ken'};
monkey_id = lower(string_fname(1:3));
monkey_id = find(ismember(monkey_names, monkey_id));
session_id = str2num(string_fname(4:6));
task_id = str2num(string_fname(8));
end