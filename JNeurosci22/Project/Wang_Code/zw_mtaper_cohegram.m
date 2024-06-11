function [c, f] = zw_mtaper_cohegram(sig1, sig2,fs, fmin, fmax)
win_size = 0.5;
bw = 8;
%% Spectrogram inspection options - Balbir Singh
params.Fs       = fs; % sampling frequency
params.tapers   = [win_size*bw, 2*win_size*bw - 2]; % taper parameters
params.pad      = -1; % pad factor for fft
params.err      = [2 0.05];
params.trialave = 0;
params.fpass = [fmin, fmax];
movingwin       = [win_size 1/fs]; %% moving window
%%
sig1 = sig1';
sig2 = sig2';
%%
[c,~,~,~,~,t,f] = cohgramc(sig1, sig2, movingwin, params);
c = permute(c, [3, 2, 1]);
c = cat(3, zeros(size(c, 1), size(c, 2), win_size*fs/2), c, zeros(size(c, 1), size(c, 2), win_size*fs/2 - 1));
end