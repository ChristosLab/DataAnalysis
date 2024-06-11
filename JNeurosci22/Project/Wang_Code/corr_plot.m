function median_corr_plot(inputs, f_range, h, labels, color_map, varargin)
yl = [-0.05, 0.3];
h_color  = [[0.2, 0.2, 0.6]; [0.6, 0.2, 0.2]; [0, 0, 0]];
h_marker = {'.', '.', '*'};
h_offset = [-.99, -1.01, -0.7];
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
    fig_size = [3, 2.5];
end

h = boolean(h);

figure('Unit', 'inches', 'Position', [2, 2, fig_size]);
hold on
for i = 1:numel(inputs)
    p(i) = plot(f_range, nanmean(inputs{i}, 1), 'DisplayName', [labels{i}, sprintf(' (N = %d)', size(inputs{i}, 1))], 'Color', color_map{i}, 'LineWidth', 1.2);
    if min(size(inputs{i})) > 1
        fill(...
            [f_range, fliplr(f_range)], ...
            [nanmedian(inputs{i}, 1) + nanstd(inputs{i}, 1)/sqrt(size(inputs{i}, 1)), fliplr(nanmean(inputs{i}, 1) - nanstd(inputs{i}, 1)/sqrt(size(inputs{i}, 1)))], color_map{i}, 'FaceAlpha', 0.15, 'EdgeAlpha', 0.15);
    end
end
for i = 1:size(h, 1)
    h_ = h(i, :);
    plot(f_range(h_), h_(h_) + h_offset(i), h_marker{i}, 'MarkerSize', 5, 'Color', h_color(i, :));
end
legend(p)
ylabel('Correlation coefficient')
xlabel('Frequency (Hz)')
xlim([f_range(1), f_range(end)]);
ylim(yl);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',font_size,'FontWeight','Bold', 'LineWidth', 1);
xlim([f_range(1), f_range(end)])
% l1 = fill([0, 0.5, 0.5, 0], [yl(1), yl(1), yl(2), yl(2)], [.9, .9, 0.2], 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
% l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
% set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% ylim([yl(1), yl(2)+0.04]);
end
    