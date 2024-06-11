function out_array = ref_nested_cell(ref, in_cell)
out_array = zeros(numel(in_cell), 1);
for i=1:numel(in_cell)
    out_array(i) = in_cell{i}{sub2ind_w(size(in_cell{i}), ref)};
end
end