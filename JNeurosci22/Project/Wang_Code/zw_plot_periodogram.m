function lines = zw_plot_periodogram(p, f, color)
l1 = mean(p, 1);
l2 = l1 + std(p, 1)./sqrt(size(p, 1));
l3 = l1 - std(p, 1)./sqrt(size(p, 1));
lines = [];
lines(1) = plot(f, l1, '-^', 'Color', color, 'MarkerSize', .5);
h = fill([f flip(f)], [l2 flip(l3)], color, 'LineStyle', 'none');
h.FaceAlpha = 0.25;
% lines(2) = plot(f, l2, '--', 'Color', color);
% lines(3) = plot(f, l3, '--', 'Color', color);
set(lines, 'LineWidth', 3);
xlabel('Passband central frequency (Hz)');
ylabel('Relative passband power');
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');
