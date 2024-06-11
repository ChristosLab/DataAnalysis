apm_files = dir(fullfile(apm_database, '*/', '*.apm'));
mat_files = dir(fullfile(apm_database, '*/', '*.mat'));
%%
task_versions = {};
for i = 1:numel(mat_files)
    load(fullfile(mat_files(i).folder, mat_files(i).name))
    task_versions = [task_versions, AllData.version];
end
%%
% clear mat_repo
for i = 1:numel(mat_files)
    load(fullfile(mat_files(i).folder, mat_files(i).name))
    mat_repo(i).parameters = AllData.parameters;
    mat_repo(i).version = AllData.version;
    mat_repo(i).Filename = lower(mat_files(i).name(1:(end - 4)));
    clear AllData
end
%%
fname_ = 'mat_file_repo';
save(fullfile(project_dir, output_database, fname_), 'mat_repo');
%%
task_durs = zeros(size(mat_repo, 2), 3);
for i = 1:size(task_durs, 1)
    try
        task_durs(i, :) = [mat_repo(i).parameters.fixationDuration, mat_repo(i).parameters.stimulusDuration, mat_repo(i).parameters.delayDuration];
    catch
        task_durs(i, :) = nan;
    end
end
%%
task_durs(strcmp({mat_repo.version}, 'ProSaccade_1dist 9-May-2012')', :)
%%  Filter out files based on task version and parameters
ps_mat_files = strcmp({mat_repo.version}, 'ProSaccade Sep-19-2011')'.*all(abs(task_durs - [1, 0.5, 1.5]) < 0.001, 2);
%%
ps_filenames = unique({mat_repo(find(ps_mat_files)).Filename});
%%
dodr_mat_files = strcmp({mat_repo.version}, 'ProSaccade_1dist 9-May-2012')';
dodr_filenames = unique({mat_repo(find(dodr_mat_files)).Filename});
%%
load('D:\Database\Zhengyang_Wang\Projects\LFP_project_v0.2\Wang_Database\lfp_tbl.mat')
%%
for i = 1:size(lfp_tbl, 1)
    lfp_tbl.ps_file(i) = any(strcmp(lfp_tbl.Filename(i), ps_filenames));
end
%%
for i = 1:size(lfp_tbl, 1)
    lfp_tbl.dodr_file(i) = any(strcmp(lfp_tbl.Filename(i), dodr_filenames));
end
%%
save('D:\Database\Zhengyang_Wang\Projects\LFP_project_v0.2\Wang_Database\lfp_tbl.mat', 'lfp_tbl')
