%%  Create dodr neuron table
% Read neuron epoch significance from Junda's spreadsheet
fname_ = fullfile(project_dir, output_database, 'NeuronEpoch_young.xlsx');
neuron1 = readtable(fname_, 'sheet', 'distractgap0102');
fname_ = fullfile(project_dir, output_database, 'NeuronEpoch_adult.xlsx');
neuron2 = readtable(fname_, 'sheet', 'distractgap0102');
fname_ = fullfile(project_dir, output_database, 'Neuron_Channels.xlsx');
neuron_channels = readtable(fname_);
% 
dodr_neuron_tbl = neuron_channels;
%
for i = 1:size(dodr_neuron_tbl, 1)
    [monkey_id_, session_id_, task_id_] = parse_fname(neuron_channels.Filename{i});
    dodr_neuron_tbl.monkey_id(i) = monkey_id_;
    dodr_neuron_tbl.session_id(i) = session_id_;
    dodr_neuron_tbl.task_id(i) = task_id_;
end
% 
new_var_range = 3:8;
for i = new_var_range
    new_var = nan(size(dodr_neuron_tbl, 1), 1);
    dodr_neuron_tbl = addvars(dodr_neuron_tbl, new_var, 'NewVariableNames', neuron1.Properties.VariableNames(i));
end
new_var = zeros(size(dodr_neuron_tbl, 1), 1);
dodr_neuron_tbl = addvars(dodr_neuron_tbl, new_var, 'NewVariableNames', 'dodr');
for neuron_ = {neuron1, neuron2}
    for i = 1:size(neuron_{1}, 1)
        [monkey_id_, session_id_, task_id_] = parse_fname(neuron_{1}.Filename{i});
        match_ = find((dodr_neuron_tbl.Neuron == neuron_{1}.Neuron(i)).*(dodr_neuron_tbl.task_id == task_id_));
        if match_ > 0
            dodr_neuron_tbl(match_, size(neuron_channels, 2) + 4:size(neuron_channels, 2) + numel(new_var_range) + 3) = neuron_{1}(i, new_var_range);
%             dodr_neuron_tbl(match_, end - numel(new_var_range) + 1:end)
%             neuron_{1}(i, new_var_range)
%             pause;
            dodr_neuron_tbl.dodr(match_) = 1;
        end
    end
end
%%
dodr_neuron_tbl = dodr_neuron_tbl(dodr_neuron_tbl.dodr == 1, :);
last_young_session = [20, 39, 91, 27];
for i = 1:size(dodr_neuron_tbl, 1)
    dodr_neuron_tbl.stage(i) = 1 + (dodr_neuron_tbl.session_id(i) > last_young_session(dodr_neuron_tbl.monkey_id(i)));
end
%%
fname_ = 'dodr_neuron_tbl.mat';
save(fullfile(project_dir, output_database, fname_), 'dodr_neuron_tbl');
%%  Create dodr neuron repo
file_names_ = strcat(dodr_neuron_tbl.Filename, '_', num2str(dodr_neuron_tbl.Neuron));
needed_neuron = ones(size(dodr_neuron_tbl, 1), 1);
if isvarname('MatData')
    clear MatData
end
for i = 1:numel(file_names_)
    try
        load(fullfile(neuron_database, file_names_{i}))
        dodr_neuron_repo(i) = MatData;
        clear MatData
    catch
        disp(['error loading neuron file: ', file_names_{i}, newline])
        needed_neuron(i) = 0;
    end
    i
end
%%
fname_ = 'dodr_neuron_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'dodr_neuron_repo');
%%

%%
function [monkey_id, session_id, task_id] = parse_fname(string_fname)
monkey_names = {'ind', 'jac', 'lem', 'ken'};
monkey_id = lower(string_fname(1:3));
monkey_id = find(ismember(monkey_names, monkey_id));
session_id = str2num(string_fname(4:6));
task_id = str2num(string_fname(8));
end