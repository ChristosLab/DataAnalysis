function tbl = zw_group_anova(temp_input)
%ZW_QUICK_ANOVA performs 2-way ANOVA on stage and task
%   Factor1: stage
%   Factor2: task
n_element = size(temp_input{1} , 2)*size(temp_input{1}, 3);
tbl = cell(6, 7, n_element);
tic
for ind = 1:n_element
    y = [temp_input{1}(:, ind); temp_input{2}(:, ind); temp_input{3}(:, ind); temp_input{4}(:, ind)]';
    y_sizes = [size(temp_input{1}, 1), size(temp_input{2}, 1), size(temp_input{3}, 1), size(temp_input{4}, 1)];
    g1 = [zeros(1, y_sizes(1)), zeros(1, y_sizes(2)), ones(1, y_sizes(3)), ones(1, y_sizes(4))];
    g2 = [ones(1, y_sizes(1)), zeros(1, y_sizes(2)), ones(1, y_sizes(3)), zeros(1, y_sizes(4))];
    [~,tbl_,~] = anovan(y, {g1, g2}, 'model', 'interaction', 'display', 'off');
    tbl(:, :, ind) = tbl_;
end
toc
sprintf('%d elements run', n_element)
end