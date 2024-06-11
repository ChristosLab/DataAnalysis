function out = max_n_ind(in_array, ord)
[~, n] = sort(in_array);
if ord < 0
    out = n(1:abs(ord));
elseif ord > 0
    out = n((end + 1 - ord):end);
end