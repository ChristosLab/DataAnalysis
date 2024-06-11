function out_array2 = zw_cwt_no_pre(...
    pre_sig, post_sig, normalizer, ...
    pre_pad, post_pad, ...
    kernel_flag, down_sample, f_range, fs...
    )
%ZW_CWT handles LFP signal segmentation, performs cwt using ZW_CMORLET_VAR
%and feeds into ZW_GROUP_CWT
%   normalizer  - 0: 1 second pre_sig
%               - 1: total post_sig
%               - 2: return rms power 
%               - 3: std (via means of MAD)
%   kernel_flag - 0: Unscaled wavelet
%               - 1: Jones wavelet
%               - 2: M.X. Cohen wavelet
%               - 3: Multi-taper from chronux

sig = [[pre_sig, post_sig]];
% Power normalization
if normalizer == 0  % Normalize amplitude by pre RMS
    sig = sig./rms(pre_sig);
elseif normalizer == 1  % Normalize amplitude by total RMS
    sig = sig./rms(sig);
elseif normalizer == 3  % Standardize signal by SD calculated from MAD
    sig = sig./(mad(sig, 1)*1.4826);
end
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
    cwt_pow2 = cwt_pow(:, (numel(pre_sig) + 1):(end - post_pad)); %    Post signal
    out_array2 = downsample_row(cwt_pow2, down_sample);
end
% pause;
end
