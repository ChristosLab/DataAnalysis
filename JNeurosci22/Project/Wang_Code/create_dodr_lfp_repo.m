%Full pipeline: Process dis-ODR task LFP files (after timestamp-aligned
%with corresponding neuron files)
%   Process files used: 
%   'mat_file_repo.mat' - created in 'mat_file_check.m'
%%
fname_ = 'mat_file_repo';
load(fullfile(project_dir, output_database, fname_), 'mat_repo');
%   Find names of dis-ODR files
%%
dodr_mat_files = strcmp({mat_repo.version}, 'ProSaccade_1dist 9-May-2012')';
dodr_filenames = unique({mat_repo(find(dodr_mat_files)).Filename});
%%
field_names = {'Filename', 'monkey_id', 'session_id', 'task_id', 'channel_id', 'stage', 'lfp_filename'};
%%
lfp_files = [dir(fullfile(lfp_database, file_identifier)); ...
    dir(fullfile(lfp_database2, file_identifier))];
monkey_id = {'ind', 'jac', 'lem', 'ken'};
last_young_session = {20, 39, 91, 27};
stage_cutoff = struct('ind', 20,'jac', 39, 'lem', 91, 'ken', 27);
%  Loops through each MAT file.
excluding_flag = 'catch';
tic
n_file = numel(lfp_files);
lfp_repo_2 = cell2struct(cell(size(field_names)), field_names, 2);
n_valid = 0;
for ind = 1:n_file
    if ~any(strcmp(lfp_files(ind).name(1:end-8), dodr_filenames))
        continue
    elseif regexp(lfp_files(ind).name, excluding_flag)
        continue
    end
    n_valid = n_valid + 1;
    lfp_repo_2(n_valid) = add_trials(lfp_repo_2(1), lfp_files(ind), stage_cutoff);
    toc
    fprintf('%d/%d\n', ind, n_file)
end
%%  Save lfp tbl
dodr_lfp_tbl = struct2table(lfp_repo_2);
fname_ = 'dodr_lfp_tbl.mat';
save(fullfile(project_dir, output_database, fname_), 'dodr_lfp_tbl');
%%  Load lfp tbl
fname_ = 'dodr_lfp_tbl.mat';
load(fullfile(project_dir, output_database, fname_));
%%  Create new lfp_repo
for i = 1:size(dodr_lfp_tbl, 1)
    if isfile(fullfile(lfp_database, dodr_lfp_tbl.lfp_filename{i}))
        load(fullfile(lfp_database, dodr_lfp_tbl.lfp_filename{i}));
    elseif isfile(fullfile(lfp_database2, dodr_lfp_tbl.lfp_filename{i}))
        load(fullfile(lfp_database2, dodr_lfp_tbl.lfp_filename{i}));
    else
        fprintf('Non-existing file %s', dodr_lfp_tbl.lfp_filename{i})
    end
    dodr_lfp_repo(i) = LFPData;
    clear LFPData
    i
end
%%
fname_ = 'dodr_lfp_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'dodr_lfp_repo');
%%  CWT parameters
fs = 500;
cue_dur = [-1 + 1/fs, 4]; % Signal length 
% sac_dur = [-3 + 1/fs, 2];
sac_dur = [-2 + 1/fs, 3]; % Edited to correct a mistake in funtion Zw_repo_cwt_new.m, see OneNote for details
normalizer = 4;
kernel_flag = 3;
down_sample = 10;
f_range = 2:2:128;
avg_method = 1; %   0 - arithmetic; 1 - geometric
%%  Compute time-frequency response
new_repo = zw_repo_cwt_new(...
    dodr_lfp_repo, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    );
%%
second_repo = zw_repo_cwt_lite(...
    dodr_lfp_repo, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    );
%%
rescue_repo = zw_repo_cwt_rescue(...
    dodr_lfp_repo, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    );
%%
rescue_counter = 0;
for i = 1:numel(new_repo)
    if numel(new_repo(i).class) == 1 && (~isempty(rescue_repo(i).class(1).norm) || numel(rescue_repo(i).class) > 1)
        new_repo(i) = rescue_repo(i);
        rescue_counter = rescue_counter + 1;
        i
    elseif numel(new_repo(i).class) == 1 && (~isempty(second_repo(i).class(1).norm) || numel(second_repo(i).class) > 1)
        new_repo(i) = second_repo(i);
        rescue_counter = rescue_counter + 1;
        i
    end
end
%%
dodr_cwt_repo = new_repo;
fname_ = 'dodr_cwt_repo';
save(fullfile(project_dir, output_database, fname_), 'dodr_cwt_repo');
%%  Compute baseline-relative time-frequency response
%%  CWT parameters
fs = 500;
cue_dur = [-1 + 1/fs, 4]; % Signal length 
sac_dur = [-2 + 1/fs, 3];
normalizer = 4;
kernel_flag = 3;
down_sample = 10;
f_range = 2:2:128;
avg_method = 1; %   0 - arithmetic; 1 - geometric
%%
cue_n_sample = (diff(cue_dur)*fs + 1)/down_sample;
sac_n_sample = (diff(sac_dur)*fs + 1)/down_sample;
%%
target_frs = [4,8; 8,16; 16, 32; 32, numel(f_range)];
baseline_bin = [26, 50];
%%
temp_cwt_cue = zeros([numel(dodr_cwt_repo), numel(f_range), cue_n_sample]);
temp_cwt_sac = zeros([numel(dodr_cwt_repo), numel(f_range), sac_n_sample]);
for i = 1:numel(dodr_cwt_repo)
    %   CueOn aligned
    if isempty(dodr_cwt_repo(i).class(1).cue_cwt)
        temp_cwt_cue(i, :, :) = nan(numel(f_range), cue_n_sample);
    else
        single_cwt_cue = zeros(0, numel(f_range), cue_n_sample);
        for j = 1:numel(dodr_cwt_repo(i).class)
            single_cwt_cue = [single_cwt_cue; dodr_cwt_repo(i).class(j).cue_cwt];
        end
        temp_cwt_cue(i, :, :) = nanmean(single_cwt_cue, 1);
    end
    %   SaccadeOn aligned
    if isempty(dodr_cwt_repo(i).class(1).sac_cwt)
        temp_cwt_sac(i, :, :) = nan(numel(f_range), sac_n_sample);
    else
        single_cwt_sac = zeros(0, numel(f_range), sac_n_sample);
        for j = 1:numel(dodr_cwt_repo(i).class)
            single_cwt_sac = [single_cwt_sac; dodr_cwt_repo(i).class(j).sac_cwt];
        end
        temp_cwt_sac(i, :, :) = nanmean(single_cwt_sac, 1);
    end
    i
end
%%
temp_cwt_cue_b = zeros([numel(dodr_cwt_repo), numel(f_range), cue_n_sample]);
temp_cwt_sac_b = zeros([numel(dodr_cwt_repo), numel(f_range), sac_n_sample]);
for i = 1:numel(dodr_cwt_repo)
    %   Cue aligned
    if isempty(dodr_cwt_repo(i).class(1).cue_cwt)
        temp_cwt_cue_b(i, :, :) = nan(numel(f_range), cue_n_sample);
    else
        single_cwt_cue = zeros(0, numel(f_range), cue_n_sample);
        for j = 1:numel(dodr_cwt_repo(i).class)
            single_cwt_cue = [single_cwt_cue; dodr_cwt_repo(i).class(j).cue_cwt./mean(dodr_cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
        end
        temp_cwt_cue_b(i, :, :) = nanmean(single_cwt_cue, 1);
    end
    %   Sac aligned
    if isempty(dodr_cwt_repo(i).class(1).sac_cwt)
        temp_cwt_sac_b(i, :, :) = nan(numel(f_range), sac_n_sample);
    else
        single_cwt_sac = zeros(0, numel(f_range), sac_n_sample);
        for j = 1:numel(dodr_cwt_repo(i).class)
            %   NOTICE: CUE aligned baseline used for SAC aligned signal
            single_cwt_sac = [single_cwt_sac; dodr_cwt_repo(i).class(j).sac_cwt./mean(dodr_cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
        end
        temp_cwt_sac_b(i, :, :) = nanmean(single_cwt_sac, 1);
    end

    i
end
%%
temp_baseline = zeros([numel(dodr_cwt_repo), numel(f_range)]);
for i = 1:numel(dodr_cwt_repo)
    if dodr_lfp_tb.task_id(i) == 2
        temp_baseline(i, :) = nan(numel(f_range), 1);
    end
    if isempty(dodr_cwt_repo(i).class(1).cue_cwt)
        temp_baseline(i, :) = nan(numel(f_range), 1);
    else
        single_cwt_cue = zeros(0, numel(f_range));
        for j = 1:numel(dodr_cwt_repo(i).class)
        single_cwt_cue = [single_cwt_cue; mean(dodr_cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
        end
        temp_baseline(i, :) = nanmean(single_cwt_cue, 1);
    end
    i
end
%%
%%  Relative cwt by class
temp_class_cwt_cue_b = zeros([numel(dodr_cwt_repo), 8, numel(f_range), cue_n_sample]);
temp_class_cwt_sac_b = zeros([numel(dodr_cwt_repo), 8, numel(f_range), sac_n_sample]);
for i = 1:numel(dodr_cwt_repo)
    %   Cue aligned
    if isempty(dodr_cwt_repo(i).class(1).cue_cwt)
        temp_class_cwt_cue_b(i, :, :, :) = nan(8, numel(f_range), cue_n_sample);
    else
        for j = 1:numel(dodr_cwt_repo(i).class)
            temp_class_cwt_cue_b(i, j, :, :) = nanmean(dodr_cwt_repo(i).class(j).cue_cwt./mean(dodr_cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3), 1);
        end
    end
    %   Sac aligned
    if isempty(dodr_cwt_repo(i).class(1).sac_cwt)
        temp_class_cwt_sac_b(i, :, :, :) = nan(8, numel(f_range), sac_n_sample);
    else
        for j = 1:numel(dodr_cwt_repo(i).class)
            %   NOTICE: CUE aligned baseline used for SAC aligned signal
            temp_class_cwt_sac_b(i, j, :, :) = nanmean(dodr_cwt_repo(i).class(j).sac_cwt./mean(dodr_cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3), 1);
        end
    end

    i
end
%%  Load manual session rejection results
fname_ = 'maunal_session_rejection_flag_dodr.mat';
load(fullfile(project_dir, output_database, fname_), 'session_flag_list', 'flag_list', 'ind_list', 'target_cwts');
%%
valid_lfp = zeros(size(dodr_lfp_tbl, 1), 1);
for i = 1:size(dodr_lfp_tbl, 1)
    if any(target_cwts(~logical(session_flag_list)) == i)
        valid_lfp(i) = 1;
    end
end   
%%%% sum(valid_lfp) = 167
%%
fname_ = 'dodr_valid_lfp.mat';
save(fullfile(project_dir, output_database, fname_), 'valid_lfp');
%%  Stage average inspection for major irregularity
% dodr_n_y = find((dodr_lfp_tbl.stage == 1).*valid_lfp);
% dodr_n_a = find((dodr_lfp_tbl.stage == 2).*valid_lfp);
dodr_n_y = find((dodr_lfp_tbl.stage == 1).*(dodr_lfp_tbl.monkey_id == 3).*valid_lfp);
dodr_n_a = find((dodr_lfp_tbl.stage == 2).*(dodr_lfp_tbl.monkey_id == 3).*valid_lfp);
%%
figure
imagesc(squeeze(nanmean(temp_cwt_cue_b(dodr_n_y,:,:), 1)) - 1, [-1,1]);
colormap(clm_)
figure
imagesc(squeeze(nanmean(temp_cwt_cue_b(dodr_n_a,:,:), 1)) - 1, [-1,1]);
colormap(clm_)
%%
figure
imagesc(squeeze(nanmean(temp_cwt_cue_b(dodr_n_y,:,:), 1))- squeeze(nanmean(temp_cwt_cue_b(dodr_n_a,:,:), 1)), [-.3,.3]);
colormap(clm_)
%%  Functions
function in_struct = add_trials(in_struct, file_dir, stage_cutoff)
%ADD_TRIALs 
%     Utility function that adds each trial in the file
%     to the trial-by-trial dataset.
%     in_struct: the trial-by-trial dataset to be written.
%     file_dir: Single file directory structure created with DIR().

monkey_id = {'ind', 'jac', 'lem', 'ken'};

data_struct_name = 'LFPData';
file_ = fullfile(file_dir.folder, file_dir.name);
assert(... %    Check for variable naming in the .MAT file.
    ismember(who('-file', file_), data_struct_name), ...
    sprintf('Data strcuture name not found: %s', data_struct_name)...
    );
%   Adds file-specific info
in_struct = add_file_info(...
    in_struct, file_dir...
    );
%   Add stage infomation - 1 young -2 adult
if in_struct.session_id > stage_cutoff.(in_struct.monkey_id)
    in_struct.stage = 2;
else
    in_struct.stage = 1;
end
in_struct.Filename = file_dir.name(1:end-8);
in_struct.lfp_filename = file_dir.name;
in_struct.monkey_id = find(ismember(monkey_id, in_struct.monkey_id));
end
%%

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