function [out_array1, out_array2] = zw_cwt(...
    pre_sig, post_sig, post_onset, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    )
%ZW_CWT handles LFP signal segmentation, performs cwt using ZW_CMORLET_VAR
%and feeds into ZW_GROUP_CWT
%   normalizer  - 0: 1 second pre_sig
%               - 1: total post_sig
%               - 2: return rms power 
%   kernel_flag - 0: Unscaled wavelet
%               - 1: Jones wavelet

% Length of data points to include after pre-fixation signal in output
pre_padding = 0.5*fs;
% Length of data points to include before post-fixation signal in output
post_padding = 0.1*fs;
% Calcuates DC component of concatanated signal
ref_mean = mean([pre_sig, post_sig]);
% Removes DC comonent from respective signals
pre_sig = pre_sig - ref_mean;
post_sig = post_sig - ref_mean;
% Calculates RMS for both 1) pre singal and 2) concatanated signal
ref_pow = [rms(pre_sig((end - 1*fs + 1):end)), rms([pre_sig, post_sig])];

% Power normalization
if normalizer == 0  % Normalize by pre RMS
    sig = [pre_sig, post_sig]./ref_pow(1);
elseif normalizer == 1  % Normalize by total RMS
    sig = [pre_sig, post_sig]./ref_pow(2);
    return
end
% Calls ZW_CMORLET_VAR to compute the average TFR of concatenated signal
cwt_amp = zw_cmorlet_var(sig, kernel_flag, f_range, fs);
% Segment the TFR back into pre and post portions for output
cwt_amp1 = cwt_amp(:, 1:(numel(pre_sig) + pre_padding)); %  Pre signal
cwt_amp2 = cwt_amp(:, ...
    (numel(pre_sig) + post_onset - 1 - post_padding):end); %    Post signal
% Downsample the TFRs by a factor of down_sample
cwt_amp1 = downsample_row(cwt_amp1, down_sample);
cwt_amp2 = downsample_row(cwt_amp2, down_sample);
% Compute the power
out_array1 = abs(cwt_amp1).^2;
out_array2 = abs(cwt_amp2).^2;

end
function outmat = downsample_row(in_mat, n)
% Quick utility function for downsampling 2D matrices by the rows without LP filtering
inmat = downsample(in_mat', n);
outmat = inmat';
end