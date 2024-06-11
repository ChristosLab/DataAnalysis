function [index, peak_time, f] = sort_trace_peak(mat_in, dur, epoch, direction, normalize, varargin)
%SORT_TRACE_PEAK finds the peaks in each row of the input matrix and sort
%the rows according to the peak timings
% Inputs:
%    mat_in    - N-by-t matrix array of time series (e.g. PSTH). with each
%                row representing one incidence and each column
%                representing one time point
%    dur       - two-element vector array specifying the total time range
%                (in unit of seconds)
%    epoch     - two-element vector array specifying the time range in which
%                to find the peak (in unit of time points)
%    direction - 0: early peaks first; 1: late peaks first
%    normalize - 1: each row is min-maxed between 0 and 1
% Outputs:
%    index     - N-by-1 vector array specifying for each original incidence
%                its new order in the sorted sequence
%    peak_time - N-by-1 vector array specifying for each original incidence
%                its peak time measured in unit of time point number
%            f - figure handle
%Example:
%   [index, peak_time, f] = sort_trace_peak(PSTH, [0, 3], [101, 200], 0, 1)
% 
%   6/2/2021 - Zhengyang Wang

%   Min-maxing each row of the input
if normalize
    mat_in = (mat_in - min(mat_in, [], 2))./(max(mat_in, [], 2) - min(mat_in, [], 2));
end
peak_time = zeros(size(mat_in, 1), 1);
%   Find peak time
for i = 1:size(mat_in, 1)
    [~, i_] = max(mat_in(i, epoch(1):epoch(2)));
    peak_time(i) = i_;

end
switch direction
    case 0
        direction = 'ascend';
    case 1
        direction = 'descend';
end
%   Sort rows by peak time
[~, index] = sort(peak_time, direction);
%   Plot
f = figure('Unit', 'inches', 'Position', [2, 2, 5, 4]);
imagesc(dur, 1:size(mat_in, 1), mat_in(index, :), [min(mat_in, [], 'all'), max(mat_in, [], 'all')]);
h = colorbar;
hold on
% ylabel(h, 'label');
end