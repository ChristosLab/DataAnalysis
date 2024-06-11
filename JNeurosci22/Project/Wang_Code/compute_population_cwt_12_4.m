%%
target_frs = [4,8; 8,16; 16, 32; 32, 64];
baseline_bin = [26, 50];
%%
temp_cwt = zeros([numel(cwt_repo), 64, 200]);
for i = 1:numel(cwt_repo)
    if lfp_tbl.task_id(i) == 2
        temp_cwt(i, :, :) = nan(64, 200);
    end
    if isempty(cwt_repo(i).class(1).cue_cwt)
        temp_cwt(i, :, :) = nan(64, 200);
    else
        single_cwt = zeros(0, 64, 200);
        for j = 1:numel(cwt_repo(i).class)
        single_cwt = [single_cwt; cwt_repo(i).class(j).cue_cwt];
        end
        temp_cwt(i, :, :) = nanmean(single_cwt, 1);
    end
    i
end
%%
save(fullfile(project_dir, output_database, 'temp_cwt.mat'), 'temp_cwt');
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt(t1, 12:17, :), 1)), 1));ylim([6, 8]*10^-3);
figure(); plot(mean(squeeze(nanmean(temp_cwt(t2, 12:17, :), 1)), 1));ylim([6, 8]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt(t1, 32:64, :), 1)), 1));ylim([.7, 1.4]*10^-3);
figure(); plot(mean(squeeze(nanmean(temp_cwt(t2, 32:64, :), 1)), 1));ylim([.7, 1.4]*10^-3);
%%
temp_cwt_b = zeros([numel(cwt_repo), 64, 200]);
for i = 1:numel(cwt_repo)
    if lfp_tbl.task_id(i) == 2
        temp_cwt_b(i, :, :) = nan(64, 200);
    end
    if isempty(cwt_repo(i).class(1).cue_cwt)
        temp_cwt_b(i, :, :) = nan(64, 200);
    else
        single_cwt = zeros(0, 64, 200);
        for j = 1:numel(cwt_repo(i).class)
        single_cwt = [single_cwt; cwt_repo(i).class(j).cue_cwt./mean(cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
        end
        temp_cwt_b(i, :, :) = nanmean(single_cwt, 1);
    end
    i
end
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_b(t1, 12:17, :), 1)), 1))%;ylim([6, 8]*10^-3);
figure(); plot(mean(squeeze(nanmean(temp_cwt_b(t2, 12:17, :), 1)), 1))%;ylim([6, 8]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_b(t1, 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
figure(); plot(mean(squeeze(nanmean(temp_cwt_b(t2, 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
%%
temp_baseline = zeros([numel(cwt_repo), 64]);
for i = 1:numel(cwt_repo)
    if lfp_tbl.task_id(i) == 2
        temp_baseline(i, :) = nan(64, 1);
    end
    if isempty(cwt_repo(i).class(1).cue_cwt)
        temp_baseline(i, :) = nan(64, 1);
    else
        single_cwt = zeros(0, 64);
        for j = 1:numel(cwt_repo(i).class)
        single_cwt = [single_cwt; mean(cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
        end
        temp_baseline(i, :) = nanmean(single_cwt, 1);
    end
    i
end
%%
fname_ = 'tfr_baseline.mat';
save(fullfile(project_dir, output_database, fname_), 'temp_baseline', 'temp_cwt_b');
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt(t1, 12:17, :)./temp_baseline(t1, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t2, 12:17, :)./temp_baseline(t2, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
figure(); plot(mean(squeeze(nanmean(temp_cwt(t1, 32:64, :)./temp_baseline(t1, 32:64), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t2, 32:64, :)./temp_baseline(t2, 32:64), 1)), 1));%ylim([6, 8]*10^-3)
figure(); plot(mean(squeeze(nanmean(temp_cwt(t1, 8:16, :)./temp_baseline(t1, 8:16), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t2, 8:16, :)./temp_baseline(t2, 8:16), 1)), 1));%ylim([6, 8]*10^-3)
%%
