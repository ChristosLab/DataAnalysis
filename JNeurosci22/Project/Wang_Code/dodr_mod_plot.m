function dodr_mod_plot(inputs, t, labels, color_map, varargin)
if numel(varargin) == 1
    yl = varargin{1};
else
    yl = [0.7, 2.3];
end
figure('Unit', 'inches', 'Position', [2, 2, 2.5, 2.5]);
hold on
for i = 1:numel(inputs)
    p(i) = plot(t, nanmean(inputs{i}, 1), 'DisplayName', [labels{i}, sprintf(' (N = %d)', size(inputs{i}, 1))], 'Color', color_map{i}, 'LineWidth', 1.5);
    fill(...
        [t, fliplr(t)], ...
        [nanmean(inputs{i}, 1) + nanstd(inputs{i}, 1)/sqrt(size(inputs{i}, 1)), fliplr(nanmean(inputs{i}, 1) - nanstd(inputs{i}, 1)/sqrt(size(inputs{i}, 1)))], color_map{i}, 'FaceAlpha', 0.15, 'EdgeAlpha', 0.15);
end
legend(p)
ylabel('Relative Power')
ylim(yl);
xlim([t(1), t(end)])
l1 = fill([0, 0.5, 0.5, 0], [yl(1), yl(1), yl(2), yl(2)], [.9, .9, 0.2], 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
l3 = fill([1, 1.5, 1.5, 1], [yl(1), yl(1), yl(2), yl(2)], [.9, .9, 0.2], 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
end
    