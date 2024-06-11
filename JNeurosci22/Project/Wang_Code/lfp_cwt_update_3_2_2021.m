%   Updates the content of lfp_repo according to the newest set of single
%   LFP files. This update is due to a potential change in time stamp
%   matching. Refrence LFP table: "lfp_tbl.mat" (The table order is
%   inherited to maintain backward compatibility). 
%   FOR FUTURE TOTAL UPDATE:
%   The current LFP files *DIFFER* from the LFP table due to the following
%   reasons:
%       1) Scope of APM file differs
%       2) Current LFP are only those with matching neuron files
%%
close all
clear
zw_setpath
fname_ = 'lfp_tbl.mat';
load(fullfile(project_dir, output_database, fname_));
%%  Create new lfp_repo
for i = 1:size(lfp_tbl, 1)
    current_filename_ = fullfile(lfp_database, lfp_tbl.lfp_filename{i});
    if isfile(current_filename_)
        load(current_filename_);
        lfp_repo(i) = LFPData;
        clear LFPData
    else
        continue
    end
    i
end
%%
clear fname_
fname_ = 'lfp_repo_3_2_2021.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_repo');
%%
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
new_repo = zw_repo_cwt_new(...
    lfp_repo, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    );
%%
cwt_repo = new_repo;
clear new_repo
clear fname_
fname_ = 'cwt_repo_3_2_2021.mat';
save(fullfile(project_dir, output_database, fname_), 'cwt_repo');
%%
rescue_repo = zw_repo_cwt_rescue(...
    lfp_repo, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    );
%%
rescue_counter = 0;
for i = 1:numel(cwt_repo)
    if numel(cwt_repo(i).class) == 1 && numel(rescue_repo(i).class) > 1
        cwt_repo(i) = rescue_repo(i);
        rescue_counter = rescue_counter + 1;
        i
    end
end
%%
clear fname_
fname_ = 'complete_cwt_repo_3_2_2-21.mat';
save(fullfile(project_dir, output_database, fname_), 'cwt_repo');