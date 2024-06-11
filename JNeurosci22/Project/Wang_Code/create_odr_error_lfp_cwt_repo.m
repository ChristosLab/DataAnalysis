%   Create lfp_error_repo from LFP files corresponding to significant
%   neurons included in the paper. Refrence LFP table: "lfp_tbl.mat" (The
%   table order is inherited to maintain backward compatibility).
% 
%   LFP files stored in directory defined by "lfp_database_error" "lfp_tbl"
%   saved in file 'lfp_tbl.mat' stored in directory defined by
%   "output_database"
%%
close all
clear
zw_setpath
fname_ = 'lfp_tbl.mat';
load(fullfile(project_dir, output_database, fname_));
%%  Create new lfp_error_repo
for i = 1:size(lfp_tbl, 1)
    [~, fn_, ext_] = fileparts(lfp_tbl.lfp_filename{i});
    %   Example error trial LFP file name: ind012_1_CH1_error.mat
    current_filename_ = fullfile(lfp_database_error, [fn_, '_error', ext_]);
    if isfile(current_filename_)
        load(current_filename_);
        if ~isempty(LFPData)
            lfp_error_repo(i) = LFPData;
        end
        clear LFPData
    else
        continue
    end
    i
end
%%
e_counter = 0;
for i = 1:numel(lfp_error_repo)
    if ~isempty(lfp_error_repo(i).class)
        e_counter = e_counter + 1;
    end
end
e_counter
% e_counter = 170
%%
clear fname_
fname_ = 'lfp_error_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_error_repo');
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
%%
cwt_error_repo = zw_repo_cwt_new(...
    lfp_error_repo, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    );
%%
% cwt_repo = new_repo;
% clear new_repo
% clear fname_
% fname_ = 'cwt_repo_3_2_2021.mat';
% save(fullfile(project_dir, output_database, fname_), 'cwt_repo');
%%
rescue_repo = zw_repo_cwt_rescue(...
    lfp_error_repo, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    );
%%
rescue_counter = 0;
for i = 1:numel(cwt_error_repo)
    if numel(cwt_error_repo(i).class) == 1 && (~isempty(rescue_repo(i).class(1).norm) || numel(rescue_repo(i).class) > 1)
        cwt_error_repo(i) = rescue_repo(i);
        rescue_counter = rescue_counter + 1;
        i
    end
end
%%
clear fname_
fname_ = 'cwt_error_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'cwt_error_repo');