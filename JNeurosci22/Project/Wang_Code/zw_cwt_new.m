function [out_array_cue, out_array_sac, norm] = zw_cwt_new(...
    pre_sig, post_sig, normalizer, ...
    cue_sample, sac_sample, ...
    kernel_flag, down_sample, f_range, fs...
    )
%ZW_CWT handles LFP signal segmentation, performs cwt using ZW_CMORLET_VAR
%and feeds into ZW_GROUP_CWT
%   normalizer  - 0: 1 second pre_sig rms
%               - 1: total post_sig rms
%               - 2: return rms power 
%               - 3: total post_sig std (via means of MAD)
%               - 4: 1 second post_sig (baseline) std (via means of MAD)
%   kernel_flag - 0: Unscaled wavelet
%               - 1: Jones wavelet
%               - 2: M.X. Cohen wavelet
%               - 3: Multi-taper from chronux
baseline_dur = [1/fs, 1]*fs;
sig = [[pre_sig, post_sig]];
baseline_sig = post_sig(baseline_dur(1):baseline_dur(2));
% Power normalization
if normalizer == 0  % Normalize amplitude by pre RMS
    norm = rms(pre_sig);
elseif normalizer == 1  % Normalize amplitude by total RMS
    norm = rms(sig);
elseif normalizer == 3  % Standardize signal by SD calculated from MAD
    norm = (mad(sig, 1)*1.4826);
elseif normalizer == 4
    norm = (mad(baseline_sig, 1)*1.4826);
end
sig = sig./norm;
if kernel_flag <3
    % Calls ZW_CMORLET_VAR to compute the average TFR of concatenated signal
    cwt_amp = zw_cmorlet_var(sig, kernel_flag, f_range, fs);
    % Segment and trim the TFR back into pre and post portions for output
    cwt_amp2 = cwt_amp(:, (numel(pre_sig) + 1):(end - post_pad)); %    Post signal
    % Downsample the TFRs by a factor of down_sample
    cwt_amp2 = downsample_row(cwt_amp2, down_sample);
    % Compute the power
    out_array2 = abs(cwt_amp2).^2;
elseif kernel_flag == 3
    cwt_pow = zw_mtaper(sig, f_range([1,end]), fs);
    cwt_pow_cue = cwt_pow(:, numel(pre_sig) + [cue_sample(1):cue_sample(2)]);
    out_array_cue = downsample_row(cwt_pow_cue, down_sample);
    cwt_pow_sac = cwt_pow(:, numel(pre_sig) + [sac_sample(1):sac_sample(2)]);
    out_array_sac = downsample_row(cwt_pow_sac, down_sample);
%     cwt_pow2 = cwt_pow(:, (numel(pre_sig) + 1):(end - post_pad)); %    Post signal
%     out_array2 = downsample_row(cwt_pow2, down_sample);
end
% pause;
end
