function out_array = zw_group_power_change(in_array, repo, interval_post, interval_pre, f_range, fs)
n_repo = numel(repo);
out_array = zeros([n_repo, size(squeeze(in_array))]);
for i = 1:numel(repo)
    if and(...
            numel(repo(i).LFP) >= diff(interval_post)*fs, ...
            numel(repo(i).pretrial_LFP) >= diff(interval_pre)*fs...
            )
        post_onset = floor((repo(i).Cue_onT + interval_post(1))*fs);
        post_offset = post_onset + diff(interval_post)*fs - 1;
        pre_onset = floor(interval_pre(1)*fs) + 1;
        pre_offset  = floor(interval_pre(2)*fs);
        post_sig = repo(i).LFP(post_onset:post_offset);
        pre_sig = repo(i).pretrial_LFP(end + pre_onset:end+ pre_offset);
        out_array(i, :, :) = zw_power_change_from_baseline(...
            pre_sig, post_sig, f_range, fs...
            );
    else
        out_array(i, :, :) = nan;
    end
end
out_array = squeeze(nanmean(out_array, 1));
end