function out = zw_tfr_patch(frequncy_band, epoch, input_tfr)
% Calculates average TFR across a time X frequency patch
% epoch - 2-dimensional vector specifying start and end sample number
% frequency_band - 2-dimensional vector specifying start and end frequency indices
if or(any(isnan(input_tfr)), numel(input_tfr) < 2)
    out = nan;
else
    out = mean(mean(input_tfr(frequncy_band(1):frequncy_band(2), epoch(1):epoch(2))));
end