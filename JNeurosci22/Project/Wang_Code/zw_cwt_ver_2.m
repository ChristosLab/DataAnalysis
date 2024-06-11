function [out_array1, out_array2] = zw_cwt_ver_2(...
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
    cwt_amp1 = cwt_amp(:, (pre_pad + 1):numel(pre_sig)); %  Pre signal
    cwt_amp2 = cwt_amp(:, (numel(pre_sig) + 1):(end - post_pad)); %    Post signal
    % Downsample the TFRs by a factor of down_sample
    cwt_amp1 = downsample_row(cwt_amp1, down_sample);
    cwt_amp2 = downsample_row(cwt_amp2, down_sample);
    % Compute the power
    out_array1 = abs(cwt_amp1).^2;
    out_array2 = abs(cwt_amp2).^2;
elseif kernel_flag == 3
    cwt_pow = zw_mtaper(sig, f_range([1,end]), fs);
    cwt_pow1 = cwt_pow(:, (pre_pad + 1):numel(pre_sig)); %  Pre signal
    cwt_pow2 = cwt_pow(:, (numel(pre_sig) + 1):(end - post_pad)); %    Post signal
    out_array1 = downsample_row(cwt_pow1, down_sample);
    out_array2 = downsample_row(cwt_pow2, down_sample);
end
% mean_pow_baseline = mean(out_array1, 2);
mean_pow_baseline = geomean(out_array1, 2);
out_array1 = out_array1./mean_pow_baseline;
% disp(mean(out_array1, 1))
% disp(mean(out_array1, 1))
% disp(mean(mean(out_array1)))
% disp(geomean(geomean(out_array1)))
out_array2 = out_array2./mean_pow_baseline;
% pause;
end
