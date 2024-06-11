zw_setpath
%%
fname_ = 'lfp_repo_complete.mat';
load(fullfile(project_dir, output_database, fname_));
%%  CWT parameters
fs = 500;
cue_dur = [-1 + 1/fs, 3]; % Signal length 
sac_dur = [-3 + 1/fs, 1];
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
%%
fname_ = 'complete_cwt_repo.mat'
save(fullfile(project_dir, output_database, fname_), 'cwt_repo');
%%
rescue_repo = zw_repo_cwt_rescue(...
    lfp_repo, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    );
%%
fname_ = 'rescue_cwt_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'rescue_repo');
%%
for i = 1:numel(cwt_repo)
    if numel(cwt_repo(i).class) == 1 && numel(rescue_repo(i).class) > 1
        cwt_repo(i) = rescue_repo(i);
        i
    end
end
%%
fname_ = 'complete_rescued_cwt_repo.mat'
save(fullfile(project_dir, output_database, fname_), 'cwt_repo');