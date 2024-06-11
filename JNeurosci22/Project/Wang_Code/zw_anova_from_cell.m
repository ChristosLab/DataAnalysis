function out = zw_anova_from_cell(anova_cell)
% 2-D unbalanced ANOVA from a N-by-M cell
y = [anova_cell{:}];
label = cell(0);
current_label_1 = cell(size(anova_cell));
current_label_2 = cell(size(anova_cell));
for j = 1:size(anova_cell, 1)
    for k = 1:size(anova_cell, 2)
        current_label_1{j, k} = zeros(size(anova_cell{j, k})) + j;
        current_label_2{j, k} = zeros(size(anova_cell{j, k})) + k;
    end
end
label{1} = [current_label_1{:}];
label{2} = [current_label_2{:}];
[~,out] = anovan(y, label, 'model', 'interaction');
end