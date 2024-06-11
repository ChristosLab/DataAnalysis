function out = zw_band_power(frequncy_band, input_tfr)
% Calculates average TFR across a time X frequency patch
% epoch - 2-dimensional vector specifying start and end sample number
% frequency_band - 2-dimensional vector specifying start and end frequency indices
if or(any(isnan(input_tfr)), numel(input_tfr) < 2)
    out = nan;
else
    out = mean(input_tfr(frequncy_band(1):frequncy_band(2), :), 1);
end