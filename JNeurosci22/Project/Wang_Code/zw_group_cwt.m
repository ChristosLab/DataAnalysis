function [out_array1, out_array2] = zw_group_cwt(...
    in_array1, in_array2, repo, ...
    interval_pre, interval_post, ...
    kernel_flag, down_sample, f_range, fs...
    )
%ZW_GROUP_CWT 1) concatenate pre- and post- cue singals; 2) calls ZW_CWT to
%compute the TFR of each trial and 3) computes the session mean of power
%density
%   in_array1       : pre cue TFR as a format template
%   in_array2       : post cue TFR as a format template
%   repo            : LFP_repo
%   interval_pre    : Interval of pre-cue signal to use in seconds, 0 =
%                     trial end
%   interval_post   : Interval of pre-cue signal to use in seconds, 0 = 
%                     cue onset
%   kernel_flag     : See ZW_CMORLET_VAR
%   down_sample     : Downsamlping factor applied to the final output in
%                     order to save memory
%   target_f_range  : Frenquencies to examine in Hertz
%   fs              : Sampling frequency in Hertz

%%
n_repo = numel(repo);
out_array1 = zeros([n_repo, size(squeeze(in_array1))]);
out_array2 = zeros([n_repo, size(squeeze(in_array2))]);

for i = 1:numel(repo)
    % Check if both the current and the previous trial are long enough
    if and(...
            numel(repo(i).LFP) >= diff(interval_post)*fs, ...
            numel(repo(i).pretrial_LFP) >= diff(interval_pre)*fs...
            )
        % Convert 
        post_onset = floor((repo(i).Cue_onT + interval_post(1))*fs);
        post_offset = post_onset + diff(interval_post)*fs - 1;
        pre_onset = floor(interval_pre(1)*fs) + 1;
        pre_offset  = floor(interval_pre(2)*fs);
        % 
        post_sig = repo(i).LFP(1:post_offset)';
        pre_sig = repo(i).pretrial_LFP(end + pre_onset:end+ pre_offset)';
        % Calls ZW_CWT for each trial
        [out_array1(i, :, :), out_array2(i, :, :)] = zw_cwt(...
            pre_sig, post_sig, post_onset, 0, ...
            kernel_flag, down_sample, f_range, fs...
            );
        
    else % Trial output assigned NaN if too short
        out_array1(i, :, :) = nan;
        out_array2(i, :, :) = nan;

    end
end
% Comuptes the session mean and trims singleton dimensions
out_array1 = squeeze(nanmean(out_array1, 1)); % 
out_array2 = squeeze(nanmean(out_array2, 1));
end