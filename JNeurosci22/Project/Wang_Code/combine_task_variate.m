function [referenced_cells_out, referenced_subs_out] = combine_task_variate(referenced_cells_in, referenced_subs_in)
% Reduce task classes to 1-8
base_class_n = 8;
reduced_subs = mod(referenced_subs_in - 1, base_class_n) + 1;
referenced_subs_out = unique(reduced_subs);
referenced_cells_out = cell(numel(referenced_subs_out), 1);
for i = 1:numel(referenced_cells_out)
    % Retain cell content size
    referenced_cells_out{i} = zeros(size(referenced_cells_in{1}));
end
for i = 1:numel(reduced_subs)
    id_ = find(referenced_subs_out == reduced_subs(i));
    referenced_cells_out{id_} = or(referenced_cells_out{id_}, referenced_cells_in{i});
end
end