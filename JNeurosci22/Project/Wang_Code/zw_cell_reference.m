function [out_cells, out_subs] = zw_cell_reference(refs, in_cells)
% Extract from an N-D cell array according to (N-1) specified dimension
% indices to output a 1-d cell array containing all non-empty cells of
% the unspecified dimension.
% Inputs:
%       - refs: N-dimensional vector with N equaling the dimensionality of
%       in_cells; 0 represents the unspecified index;
%       - groups: N-D cell array
% Outputs:
%       - out_cells: Non-empty cells from the groups
%       - out_subs: Reference subscripts of the corresponding non-empty
%       cells
cell_dims = size(in_cells);
subs = refs + (refs == 0).*(1:cell_dims(refs == 0))';
out_cells = cell(size(subs, 1), 1);
out_label_inds = 1:size(subs, 1);
nonempty_filter = [];
for i = out_label_inds
    out_cells{i} = in_cells{sub2ind_w(cell_dims, subs(i, :))};
    if sum(out_cells{i}) ~= 0 % Assuming the content of each cell is a boolean array
        nonempty_filter = [nonempty_filter, i];
    end
end
out_cells = out_cells(nonempty_filter);
if numel(out_cells) == 0
    out_subs = [];
else
    out_subs = out_label_inds(nonempty_filter);
%     out_label_inds = out_label_inds(nonempty_filter);
%     out_subs = refs + (refs == 0).*out_label_inds';
end
end