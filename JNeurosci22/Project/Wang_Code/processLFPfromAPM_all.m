clear
close all
zw_setpath;
%%  Define task id
required_task_id = [1, 2, 3 ,4 , 5];
%%  Load neurons included in the study
load('D:\Database\Zhengyang_Wang\Projects\LFP_project_v0.2\Wang_Database\neuron_tbl.mat');
%%  Find out already prsent files
lfp_files = dir(fullfile(lfp_database, file_identifier));
present_files = cell(numel(lfp_files), 1);
for i = 1:numel(present_files)
    [~, present_files{i}, ~] = fileparts(lower(lfp_files(i).name));
    present_files{i} = present_files{i}(1:8);
end
present_files = unique(present_files);
%%  Step1: Process Prosaccade files by taking timestamps from behavior file
progress_counter = 0;
processed_counter = 0;
processed_files = [];
neuron_tbl.Filename = lower(neuron_tbl.Filename);
expected_files = neuron_tbl.Filename(ismember([neuron_tbl.task_id], required_task_id));
for i = 1:numel(apm_files)
    [~, fname, ~] = fileparts(apm_files(i).name);
    fname = lower(fname);
    if any(strcmp(present_files, fname))
        disp([fname, ' already processed: skipping', newline]);
        continue
    end
    if any(strcmp(expected_files, fname))
        channels = unique(neuron_tbl.Channel(find(strcmp(neuron_tbl.Filename, fname))));
        try
            AllData = [];
            load(fullfile(apm_files(i).folder, [upper(fname),'.mat']))  %%Loads AllData
%             if ~check_task(AllData)
%                 disp(['Not the required task: ', fname])
%                 continue
%             end
        catch
            disp(['problem getting behavioral data, exiting ', fname])
            continue
        end
        for j = 1:numel(channels)
            LFPData = LFP_Data_Extraction(apm_files(i).folder, fname, channels(j));
            save(fullfile(lfp_database, [fname, '_CH', num2str(channels(j)), '.mat']), 'LFPData')
        end
        processed_counter = processed_counter + 1;
        processed_files = [processed_files, newline, fname];
        fid = fopen(fullfile(lfp_database, 'processed_files.txt'), 'wb+');
        fwrite(fid, processed_files);
        fclose(fid);
    else
        disp(['Not a neuron-associated file: ', fname])
    end
    progress_counter = i
    fid = fopen(fullfile(lfp_database, 'progress_counter.txt'), 'wb+');
    fwrite(fid, num2str(progress_counter));
    fclose(fid);
end
%%  Step2: Match timestamps from neuron data
load('D:\Database\Zhengyang_Wang\Projects\LFP_project_v0.2\Wang_Database\neuron_tbl_w_pfc.mat');
%%
lfp_files = dir(fullfile(lfp_database, file_identifier));
err_log = [];
err_mat = [];
converted_file_id = [];
for i = 1:numel(lfp_files)
    i
    %     fn_ = lfp_files(i).name;
    %     fn_(end - 5) = '';
    %     movefile(fullfile(lfp_files(i).folder, lfp_files(i).name), fullfile(lfp_files(i).folder, fn_));
    load(fullfile(lfp_files(i).folder, lfp_files(i).name));
    [field_names, field_values] = read_lfp_fname(lfp_files(i).name);
%     if strcmp(field_values{3}, '2')
%         continue
%     end
    match_ = ones(size(neuron_tbl, 1), 1);
    for j = 1:numel(field_names)
%         value_ = field_values{j};
%         if ischar(value_) && all(ismember(value_, '1234567890'))
%             field_values{j} = str2double(value_);
%         end
        match_ = match_.*strcmp(string([neuron_tbl.(field_names{j})]), string(field_values{j}));
    end
    match_ = find(match_);
    if ~isempty(match_)
        run_flag = 1;
        neuron_counter = 0;
        while run_flag && (neuron_counter < numel(match_))
            neuron_counter = neuron_counter + 1;
            if isempty(neuron_repo(match_(neuron_counter)).class)
                err_msg = ['Non-existing neuron data for neuron no.', num2str(neuron_tbl.Neuron(match_(neuron_counter))), newline];
                disp(err_msg)
                err_code = 4;
                err_log = [err_log, err_msg];
                err_mat = [err_mat; [i, match_(neuron_counter), err_code]];
                fid = fopen(fullfile(lfp_database, 'match_timestamp_err_log_total.txt'), 'wb+');
                fwrite(fid, err_log);
                fclose(fid);
                continue
            end
            [LFPData, err_code] = match_timestamp(LFPData, neuron_repo(match_(neuron_counter)));
            switch err_code
                case {0, 3}
                    save(fullfile(lfp_files(i).folder, lfp_files(i).name), 'LFPData');
                    run_flag = 0;
                    err_msg = [];
                    converted_file_id =  [converted_file_id, i];
                    disp([num2str(numel(converted_file_id)), ' of files  converted', newline])
                    if err_code == 3
                        err_msg = ['Some trial numbers did not match between ', lfp_files(i).name, ' and neuron no.', num2str(neuron_tbl.Neuron(match_(neuron_counter))), newline];
                    end
                case 1  %   Non-matching trial numbers
                    err_msg = ['No matching trial numbers between LFP file: ', lfp_files(i).name, ' and neuron no.', num2str(neuron_tbl.Neuron(match_(neuron_counter))), newline];
                case 2
                    err_msg = ['Cue on time too different between LFP file: ', lfp_files(i).name, ' and neuron no.', num2str(neuron_tbl.Neuron(match_(neuron_counter))), newline];
            end
            if ~isempty(err_msg)
                disp(err_msg)
                err_log = [err_log, err_msg];
                err_mat = [err_mat; [i, match_(neuron_counter), err_code]];
                fid = fopen(fullfile(lfp_database, 'match_timestamp_err_log_total.txt'), 'wb+');
                fwrite(fid, err_log);
                fclose(fid);
            end
            
        end
    else
        err_msg = ['No neurons matched for ', lfp_files(i).name, newline];
        disp(err_msg);
        err_log = [err_log, err_msg];
        err_mat = [err_mat; [i, 0, err_code]];
        fid = fopen(fullfile(lfp_database, 'match_timestamp_err_log_total.txt'), 'wb+');
        fwrite(fid, err_log);
        fclose(fid);
    end
end
%%
save(fullfile(lfp_files(i).folder, 'err_mat_total_updated.mat'), 'err_mat');
%%
check_task(AllData_681)
%%
function out = check_task(mat_data)
out = (abs(mat_data.parameters.fixationDuration - 1) < 10^-6)* ...
    (abs(mat_data.parameters.stimulusDuration - 0.5) < 10^-6)* ...
    (abs(mat_data.parameters.delayDuration -1.5) < 10^-6);
end
%%
function [field_names, field_values] = read_lfp_fname(fname)
%READ_LFP_FNAME
%     Utility function that deliminate file name strings in LFP saving
%     conventions. Feeds into ADD_FILE_INFO

field_names = {'monkey_id', 'session_id', 'task_id', 'Channel'};
field_values = textscan(fname, '%s%s%s', 'Delimiter', '_'); % Deliminate FNAME using "_".
ch_scan = textscan(field_values{3}{1}, '%2c%f');
field_values = [textscan(field_values{1}{1}, '%3c%f'), field_values{2}{1}, ch_scan{2}]; % Seperate SUBJECT_ID and SESSION_ID.
% Convert monkey id to dummy code
monkey_id = {'ind', 'jac', 'lem', 'ken'};
field_values{1} = find(strcmp(monkey_id, field_values(1)));
end