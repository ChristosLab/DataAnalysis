function [out_array1, out_array2] = zw_group_cwt_ver_2(...
    in_array1, in_array2, repo, ...
    pre_dur, post_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs, ...
    avg_method...
    )
%ZW_GROUP_CWT_VER_2 1) concatenate pre- and post- cue singals; 2) calls ZW_CWT to
%compute the TFR of each trial and 3) computes the session mean of power
%density
%  ver_2 created on 6_2_20: Pre-fixation period is now decided after the
%  pre- and current trial are concatenated; courtesy of Balbir Singh.

%   in_array1       : pre cue TFR as a format template
%   in_array2       : post cue TFR as a format template
%   repo            : LFP_repo
%   pre_dur         : Duration before cue onset
%   post_dur        : Duration starting from cue onset
%   kernel_flag     : See ZW_CMORLET_VAR
%   down_sample     : Downsamlping factor applied to the final output in
%                     order to save memory
%   target_f_range  : Frenquencies to examine in Hertz
%   fs              : Sampling frequency in Hertz
%   avg_method   - 0: Arithmetic mean
%                - 1: Geometric mean


%%  Pad to avoid edge effects
pre_pad = 0.3 * fs; %  number of samples to pad towards pre-trial
post_pad = 1 * fs; % number pf samples length to pad towards post-trial
%%
n_repo = numel(repo);
out_array1 = zeros([n_repo, size(squeeze(in_array1))]);
out_array2 = zeros([n_repo, size(squeeze(in_array2))]);

parfor i = 1:numel(repo)
% parfor i = 1:numel(repo)
    % Check if both the current and the previous trial are long enough
    pre_trial_T = repo(i).Cue_onT - 1;
    pre_sample_ava = numel(repo(i).pretrial_LFP) + floor(pre_trial_T * fs);
    pre_sample_req = pre_pad + pre_dur*fs; % Pad backward in time
    post_sample_ava = numel(repo(i).LFP) - floor(pre_trial_T * fs) ;
    post_sample_req = post_pad + post_dur*fs; %  Pad forward in time
    if and(...
            pre_sample_ava > pre_sample_req, ...
            post_sample_ava > post_sample_req...
            )
        % Convert 
        sig = [repo(i).pretrial_LFP', repo(i).LFP'];
        % Filtering
        sig = zw_lfp_filt(sig);
        pre_sig = sig((pre_sample_ava - pre_sample_req + 1):(pre_sample_ava));
        post_sig = sig((pre_sample_ava + 1):(pre_sample_ava + post_sample_req));
        % Calls ZW_CWT_VER_2 for each trial
        [o1, o2] = zw_cwt_ver_2(...
            pre_sig, post_sig, normalizer, ...
            pre_pad, post_pad, ...
            kernel_flag, down_sample, f_range, fs...
            );
        out_array1(i, :, :) = o1;
        out_array2(i, :, :) = o2;
    else % Trial output assigned NaN if too short
        out_array1(i, :, :) = nan;
        out_array2(i, :, :) = nan;

    end
end
% Comuptes the session mean and trims singleton dimensions
if avg_method == 0
    out_array1 = squeeze(nanmean(out_array1, 1));
    out_array2 = squeeze(nanmean(out_array2, 1));
elseif avg_method == 1
    out_array1 = squeeze(geomean(out_array1, 1, 'omitnan'));
    out_array2 = squeeze(geomean(out_array2, 1, 'omitnan'));
end
end