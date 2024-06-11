%Convert ERROR TRIALS from ODR apm files included in study.
clear
close all
zw_setpath;
%%  Load LFP sites included in the study
load(fullfile(project_dir, output_database, 'lfp_tbl.mat'));
load(fullfile(project_dir, output_database, 'neuron_site_categories.mat'));
load(fullfile(project_dir, output_database, 'neuron_site_categories_nonsig.mat'));
%%  Define matching APM file names and channel for desired LFP files
target_lfp_tbl = lfp_tbl([site_mod_cat{:, :, 1:2}], :);
target_files = unique(target_lfp_tbl.Filename);
target_channels = cell(size(target_files));
for i = 1:numel(target_files)
    target_channels{i} = target_lfp_tbl.channel_id(ismember(target_lfp_tbl.Filename, target_files{i}));
end
%%  Find out already prsent error files
lfp_files = dir(fullfile(lfp_database_error, file_identifier));
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
for i = 1:numel(apm_files)
    [~, fname, ~] = fileparts(apm_files(i).name);
    fname = lower(fname);
    if any(strcmp(present_files, fname))
        disp([fname, ' already processed: skipping', newline]);
        continue
    end
    if any(strcmp(target_files, fname))
        channels = target_channels{ismember(target_files, fname)};
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
            %   12/12/21 -ZW
            %   Updated extraction function to correctly catch all saccade
            %   errror trials, previously conly capturing those having gone
            %   into the target window
%             LFPData = LFP_Data_Extraction_error(apm_files(i).folder, fname, channels(j));
            LFPData = LFP_Data_Extraction_error_edited(apm_files(i).folder, fname, channels(j));
            save(fullfile(lfp_database_error, [fname, '_CH', num2str(channels(j)), '_error.mat']), 'LFPData')
        end
        processed_counter = processed_counter + 1;
        processed_files = [processed_files, newline, fname];
        fid = fopen(fullfile(lfp_database_error, 'processed_files.txt'), 'wb+');
        fwrite(fid, processed_files);
        fclose(fid);
    else
        disp(['Not a target file: ', fname])
    end
    progress_counter = i
    fid = fopen(fullfile(lfp_database_error, 'progress_counter.txt'), 'wb+');
    fwrite(fid, num2str(progress_counter));
    fclose(fid);
end

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