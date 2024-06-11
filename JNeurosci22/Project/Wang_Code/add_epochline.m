function add_epochline(line_color)
yl = ylim;
l1 = line([0, 0], yl, 'LineStyle', '--', 'LineWidth',1.2, 'Color', line_color);
l2 = line([.5, .5], yl, 'LineStyle', '--', 'LineWidth',1.2, 'Color', line_color);
l3 = line([2, 2], yl, 'LineStyle', '--', 'LineWidth',1.2, 'Color', line_color);
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
end
