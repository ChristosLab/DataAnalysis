function [out, stats] = zw_anova1_from_cell(anova_cell)
% 1-D unbalanced ANOVA from a N-by-1 cell
y = [anova_cell{:}];
label = cell(0);
current_label_1 = cell(size(anova_cell));
for j = 1:numel(anova_cell)
        current_label_1{j} = zeros(size(anova_cell{j})) + j;
end
label{1} = [current_label_1{:}];
[~,out, stats] = anovan(y, label, 'model', 'interaction');
end