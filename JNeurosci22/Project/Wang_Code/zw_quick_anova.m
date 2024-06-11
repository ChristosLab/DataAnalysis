function tbl = zw_quick_anova(periodos_young_ps, periodos_young_as, periodos_adult_ps, periodos_adult_as, ind)
%ZW_QUICK_ANOVA performs 2-way ANOVA on stage and task
%   Factor1: stage
%   Factor2: task
y = [periodos_young_ps(:, ind);periodos_young_as(:, ind); periodos_adult_ps(:, ind); periodos_adult_as(:, ind)]';
y_sizes = [size(periodos_young_ps, 1), size(periodos_young_as, 1), size(periodos_adult_ps, 1), size(periodos_adult_as, 1)];
g1 = [zeros(1, y_sizes(1)), zeros(1, y_sizes(2)), ones(1, y_sizes(3)), ones(1, y_sizes(4))];
g2 = [ones(1, y_sizes(1)), zeros(1, y_sizes(2)), ones(1, y_sizes(3)), zeros(1, y_sizes(4))];
[~,tbl,stats] = anovan(y, {g1, g2}, 'model', 'interaction', 'display', 'off');
end
