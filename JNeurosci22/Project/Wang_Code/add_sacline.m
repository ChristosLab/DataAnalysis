function add_sacline(line_color)
yl = ylim;
l1 = line([0, 0], yl, 'LineStyle', '--', 'LineWidth',2, 'Color', line_color);
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
end
