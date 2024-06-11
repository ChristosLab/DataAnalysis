function pev_plot(inputs, t, labels, color_map, varargin)
yl = [-0.01, 0.145];
if numel(varargin) > 0
    if ~isempty(varargin{1})
        yl = varargin{1};
    end
end
if numel(varargin) > 1
    font_size = varargin{2};
else
    font_size = 10;
end
if numel(varargin) > 2
    fig_size = varargin{3};
else
    fig_size = [2.5, 2.5];
end
figure('Unit', 'inches', 'Position', [2, 2, fig_size]);
hold on
for i = 1:numel(inputs)
    p(i) = plot(t, nanmean(inputs{i}, 1), 'DisplayName', [labels{i}, sprintf(' (N = %d)', size(inputs{i}, 1))], 'Color', color_map{i}, 'LineWidth', 1.2);
    fill(...
        [t, fliplr(t)], ...
        [nanmean(inputs{i}, 1) + nanstd(inputs{i}, 1)/sqrt(size(inputs{i}, 1)), fliplr(nanmean(inputs{i}, 1) - nanstd(inputs{i}, 1)/sqrt(size(inputs{i}, 1)))], color_map{i}, 'FaceAlpha', 0.15, 'EdgeAlpha', 0.15);
end
legend(p)
ylabel('PEV')
ylim(yl);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',font_size,'FontWeight','Bold', 'LineWidth', 1);
xlim([t(1), t(end)])
l1 = fill([0, 0.5, 0.5, 0], [yl(1), yl(1), yl(2), yl(2)], [.9, .9, 0.2], 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% ylim([yl(1), yl(2)+0.04]);
end
    