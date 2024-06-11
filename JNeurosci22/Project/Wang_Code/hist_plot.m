function hist_plot(inputs, bin_edges, labels, color_map, norm_method)
n_inputs = numel(inputs);
to_plot = zeros(n_inputs, numel(bin_edges) - 1);
for i = 1:n_inputs
    to_plot(i, :) = histcounts(inputs{i}, bin_edges, 'Normalization', norm_method);
end
% 
figure('Unit', 'inches', 'Position', [2, 2, 3, 2]);
b = bar([bin_edges(1:end-1)] + diff(bin_edges)/2, to_plot', 'BarWidth', 1);
for i = 1:numel(b)
    b(i).FaceColor = color_map{i};
    b(i).DisplayName = labels{i};
end
ylabel('Probability');
xlabel('Correlation coefficient');
ylim([0, 0.18]);
legend
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
end