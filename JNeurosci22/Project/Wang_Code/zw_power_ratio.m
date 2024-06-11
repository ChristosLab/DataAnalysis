function p = zw_power_ratio(sig, lb, ub, f_interval, fs)
%ZW_POWER_RATIO Computes average power ratio of each passband 
% relative to the given range.
% Detrend
sig = sig - mean(sig);
% Normalize to total power of 1
sig = sig./rms(sig);
[p_, f_] = periodogram(sig, [], length(sig), fs);
origingal_f_interval = mean(diff(f_));
f_band = [...
    [lb:f_interval:ub - f_interval + origingal_f_interval]', ...
    [lb + f_interval - origingal_f_interval:f_interval:ub]'...
    ];
total_band = [lb, ub];
p = zeros(size(f_band, 1), 1);
for i = 1:size(f_band, 1)
    p(i) = bandpower(p_, f_, f_band(i, :), 'psd')/bandpower(p_, f_, total_band, 'psd');
% Convert power density to power.
end
end