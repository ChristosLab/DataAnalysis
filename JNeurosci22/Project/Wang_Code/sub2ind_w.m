function ind1d = sub2ind_w(dims, array_sub)
array_sub = array_sub - 1;  % Converting to zero-based indexing
ind1d = 0;
for digi = 1:numel(dims)
    ind1d = ind1d + array_sub(digi)*prod(dims(1:(digi - 1)));
end
ind1d = ind1d + 1;  % Reverting back to one-based indexing
end