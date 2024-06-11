close all;
clear;
zw_setpath;
%%  Preallocate structure field names.
field_names = {...
    'monkey_id', 'session_id', 'task_id', 'channel_id', ...
    'Cue_onT', 'Syncs', 'Trial_Num', 'Sampling_Rate', 'LFP', 'pretrial_LFP', ...
    'voltage_calibration', 'LSG', 'class', 'stage', 'pretrial_flag' ...
    };
%%
monkey_id = {'ind', 'jac', 'lem', 'ken'};
last_young_session = {20, 26, 66, 27};
stage_cutoff = struct('ind', 20,'jac', 26, 'lem', 66, 'ken', 27);
lfp_repo_2 = cell2struct(cell(size(field_names)), field_names, 2);
%  Loops through each MAT file.
excluding_flag = 'catch';
tic
n_file = numel(lfp_files);
for ind = 1:n_file
    lfp_repo_2 = add_trials(lfp_repo_2, lfp_files(ind), stage_cutoff, excluding_flag);
    toc
    sprintf('%d/%d\n', ind, n_file)
end
%%  Save 2nd batch LFP repo
save(fullfile(project_dir, output_database, 'lfp_repo_2.mat'), 'lfp_repo_2');
%%  Functions
function out_struct = add_trials(in_struct, file_dir, stage_cutoff, excluding_flag)
%ADD_TRIALs 
%     Utility function that adds each trial in the file
%     to the trial-by-trial dataset.
%     in_struct: the trial-by-trial dataset to be written.
%     file_dir: Single file directory structure created with DIR().

monkey_id = {'ind', 'jac', 'lem', 'ken'};
out_struct = in_struct;
if regexp(file_dir.name, excluding_flag)
    return
end
data_struct_name = 'LFPData';
n_entry = numel(in_struct); %   Input entry number.
file_ = fullfile(file_dir.folder, file_dir.name);
assert(... %    Check for variable naming in the .MAT file.
    ismember(who('-file', file_), data_struct_name), ...
    sprintf('Data strcuture name not found: %s', data_struct_name)...
    );
loaded_ = getfield(load(file_), data_struct_name);
%   Check pretrial availability
trial_Num_pool = [];
for i =1:numel(loaded_.class)
    if ~isempty(loaded_.class(i).ntr)
        trial_Num_pool = [trial_Num_pool, [loaded_.class(i).ntr.Trial_Num]];
    end
end
[pre, post] = zw_find_trials_with_pre(trial_Num_pool);
pretrial_LFP = cell(numel(pre), 1);
for i = 1:numel(loaded_.class)
    if ~isempty(loaded_.class(i).ntr)
        for j = 1:numel(loaded_.class(i).ntr)
            pre_check = find(pre == loaded_.class(i).ntr(j).Trial_Num);
            if pre_check
            % Check for empty recordings.
                if isempty(loaded_.class(i).ntr(j).LFP)
                    pretrial_LFP{pre_check} = [];
                else
                    pretrial_LFP{pre_check} = loaded_.class(i).ntr(j).LFP;
                end
            end
        end
    end
end

%   Iterates through classes then trials.
for i = 1:numel(loaded_.class)
    if ~isempty(loaded_.class(i).ntr)
        for j = 1:numel(loaded_.class(i).ntr)
            % Check for empty recordings.
            if isempty(loaded_.class(i).ntr(j).LFP)
                continue
            end
            %   Iterate through fields.
            for field_name_ = fieldnames(loaded_.class(i).ntr(j))'
                field_name_ = field_name_{1};
                %   Adds trial-specific info
                value_ = loaded_.class(i).ntr(j).(field_name_);                    
                if ischar(value_)
                    if all(ismember(value_, '1234567890'))
                        value_ = str2double(value_);
                    end
                end
                out_struct(n_entry).(field_name_) = value_;
            end
            %   Adds file-specific info            
            out_struct(n_entry) = add_file_info(...
                out_struct(n_entry), file_dir...
                );
            %   Add class ID in strcuture
            out_struct(n_entry).class = i;
            %   Add stage infomation - 1 young -2 adult
            if out_struct(n_entry).session_id > stage_cutoff.(out_struct(n_entry).monkey_id)
                out_struct(n_entry).stage = 2;
            else
                out_struct(n_entry).stage = 1;
            end
            out_struct(n_entry).monkey_id = find(ismember(monkey_id, out_struct(n_entry).monkey_id));
            %   Add pre-trial LFP if applicable
            pretrial_check = find(post == out_struct(n_entry).Trial_Num);
            out_struct(n_entry).pretrial_flag = 0;
            if pretrial_check
                pre_LFP_temp_ = pretrial_LFP{pretrial_check};
                out_struct(n_entry).pretrial_LFP = pre_LFP_temp_;
                if ~isempty(pre_LFP_temp_)
                    out_struct(n_entry).pretrial_flag = 1;
                end
            end
            %   Next loop
            n_entry = n_entry + 1;
        end
    end
end
end
%%
function in_struct = add_file_info(in_struct, file_dir)
%ADD_FILE_INFO 
%     Utility function that adds information in the file name
%     to the trial-by-trial dataset.
%     in_struct: the trial-by-trial dataset to be written.
%     file_dir: Single file directory structure created with DIR(). Feeds
%     into ADD_TRIAL_INFO

[field_names, field_values] = read_lfp_fname(file_dir.name); % 'name' field referenced.
% out_struct = in_struct;
for j = 1:numel(field_names)
    value_ = field_values{j};
    if ischar(value_)
        if all(ismember(value_, '1234567890'))
            value_ = str2double(value_);
        end
    end
    in_struct.(field_names{j})= value_;
end       
end
%%
function [field_names, field_values] = read_lfp_fname(fname)
%READ_LFP_FNAME
%     Utility function that deliminate file name strings in LFP saving 
%     conventions. Feeds into ADD_FILE_INFO
    
field_names = {'monkey_id', 'session_id', 'task_id', 'channel_id'};
field_values = textscan(fname, '%s%s%s', 'Delimiter', '_'); % Deliminate FNAME using "_".
ch_scan = textscan(field_values{3}{1}, '%2c%f');
field_values = [textscan(field_values{1}{1}, '%3c%f'), field_values{2}{1}, ch_scan{2}]; % Seperate SUBJECT_ID and SESSION_ID.
end