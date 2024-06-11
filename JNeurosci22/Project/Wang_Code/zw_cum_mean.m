function [percentile, means] = zw_cum_mean(input)
%PERCENTILE_INDICES computes the means of the first x% elements
%[percentile, means] = zw_cum_mean(input)
n_input = numel(input);
means = zeros(1, n_input);
percentile = (1:n_input)/n_input;
for i = 1:n_input
    means(i) = mean(input(1:i));
end
end