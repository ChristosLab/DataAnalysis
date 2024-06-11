function [c, f] = zw_mtaper_coherence(sig1, sig2, fs, fmin, fmax)
bw = 8;
fr_res = 2;
%% Spectrogram inspection options - Balbir Singh
params.Fs       = fs; % sampling frequency
params.tapers   = [size(sig1, 2)/fs, bw, size(sig1, 2)/fs*2*bw - 6]; % taper parameters
% params.tapers   = [win_size*bw, 2*win_size*bw - 2]; % taper parameters
params.pad      = -1; % pad factor for fft
params.err      = [2 0.05];
params.trialave = 0;
params.fpass = [fmin, fmax];
%%
sig1 = sig1';
sig2 = sig2';
%%
[c,~,~,~,~,f] = coherencycpb(sig1, sig2, params);
if mean(diff(f)) < fr_res
    c = downsample(c, round(fr_res/mean(diff(f))));
end
c = permute(c, [2, 1]);
end