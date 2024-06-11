function out_array = zw_stft(...
    pre_sig, post_sig, normalizer, ...
    window, overlap, f_range, fs...
    )
%ZW_STFT performs stft and feeds into ZW_GROUP_STFT
%   normalizer - 0: 1 second pre_sig
%              - 1: total post_sig
%              - 2: return rms power 
% Calcuates DC component of concatanated signal
ref_mean = mean([pre_sig, post_sig]);
% Removes DC component from respective signals
pre_sig = pre_sig - ref_mean;
post_sig = post_sig - ref_mean;
ref_pow = [rms(pre_sig((end - 1*fs + 1):end)), rms([pre_sig, post_sig])];

% Power normalization
if normalizer == 0
    sig = [pre_sig, post_sig]./ref_pow(1);
elseif normalizer == 1
    sig = [pre_sig, post_sig]./ref_pow(2);
elseif normalizer == 2
    out_array = ref_pow;
    return
end
stft = spectrogram(sig, window, overlap, f_range, fs);
out_array = abs(stft).^2;
end