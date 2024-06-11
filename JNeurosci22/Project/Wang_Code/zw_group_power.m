function p = zw_group_power(periodo, repo, dur, lb, ub, f_interval, fs)
n_repo = numel(repo);
p = zeros(n_repo, numel(periodo));
for i = 1:numel(repo)
    if numel(repo(i).LFP) >= dur*fs
        p(i, :) = zw_power_ratio(repo(i).LFP(1:dur*fs), lb, ub, f_interval, fs);
    else
        p(i, :) = nan;
    end
end
p = nanmean(p, 1);
end