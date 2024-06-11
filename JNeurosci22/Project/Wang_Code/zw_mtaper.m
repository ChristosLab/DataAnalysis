function pow = zw_mtaper(sig, fpass, fs)
win_size = 250/fs;
bw = 8;
%% Spectrogram inspection options - Balbir Singh
params.Fs       = fs; % sampling frequency
params.fpass    = fpass; % band of frequencies to be kept
params.tapers   = [win_size*bw, 2*win_size*bw - 2]; % taper parameters
params.pad      = -1; % pad factor for fft
params.err      = [2 0.05];
params.trialave = 0;
movingwin       = [win_size 1/fs]; %% moving window
%%
[pow,~,~,~] = mtspecgramc(sig, movingwin, params);
pow = pow';
pow = [zeros(size(pow, 1), win_size*fs/2), pow, zeros(size(pow, 1), win_size*fs/2 - 1)];
end