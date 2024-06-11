function out = zw_extract_entry(repo, targets, values, varargin)
%ZW_EXTRACT_ENTRY Extract structure indices with certain field values.
%     REPO must be a structure.
%     TARGETS and VALUES must be cells.
%     Example:
%     temp_s = struct();
%     temp_s(1).a = 1;
%     temp_s(2).a = 1;
%     temp_s(3).a = 1;
%     temp_s(1).b = 4;
%     temp_s(2).b = 3;
%     temp_s(3).b = 3;
% 
%     targets = {'a', 'b'};
%     values = {1,3};
% 
%     out = zw_extract_entry(temp_s, targets, values)
%     temp_s(out)
% 
%     Outputs:
%     ut =
% 
%          2
%          3
% 
% 
%     ans = 
% 
%       1×2 struct array with fields:
% 
%         a
%         b

uo_flag = true;
if numel(varargin) > 0
    uo_flag = varargin{1};
end

indices = zeros(numel(repo), numel(targets)); % Temporary indices for each condition.
for ind = 1:length(targets)
    fun = @(x) getfield(repo(x), targets{ind}) == values{ind};
    indices(:, ind) = arrayfun(fun, 1:numel(repo), 'UniformOutput', uo_flag);
end
out = find(prod(indices, 2));
end
