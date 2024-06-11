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
            check_error_assumption(AllData);
        catch
            disp(['problem getting behavioral data, exiting ', fname])
            continue
        end
    end
%     i
end
function check_error_assumption(AllData)
for i = 1:numel(AllData.trials)
%     AllData.trials(i).EndofTrialtime - AllData.trials(i).CueOff
    if (AllData.trials(i).EndofTrialtime - AllData.trials(i).CueOff) < 1.5 && AllData.trials(i).Statecode == 5
%     if ~isempty(AllData.trials(i).CueOff) && AllData.trials(i).Statecode == 4
%     if strcmp(AllData.trials(i).Reward, 'No') && AllData.trials(i).Statecode == 6
%     if isempty(AllData.trials(i).TargetIn) && (ismember(AllData.trials(i).Statecode, [5,6]))
        i
    end
end
end