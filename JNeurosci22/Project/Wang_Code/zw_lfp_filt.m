function out = zw_lfp_filt(in)
%%  Filtering options. - Balbir Singh
N       = 5; %  2Nth order bandpass filter (See document for BUTTER()); FILTFILT() further results in 4Nth order
flp     = 1;
fhi     = 200;
fs      = 500;
wn      = [flp/(fs/2), fhi/(fs/2)];
[bb,aa] = butter(N,wn,'bandpass');
pl_freq = 60;
notch_w = 1/10;
stop_band = (pl_freq - 1 * (notch_w)):(notch_w):(pl_freq + 1 * (notch_w)); % Powerline noise removal.
in = in - mean(in, 2); %   Zero-mean
filtered_ = filtfilt(bb,aa,in');
filtered_ = filtered_';
out = ft_preproc_dftfilter(filtered_,fs, stop_band);  %% power line noise removal
end