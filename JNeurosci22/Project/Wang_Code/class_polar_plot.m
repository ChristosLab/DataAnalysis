function [p_y_, p_a_] = class_polar_plot(relative_class_counts_y, relative_class_counts_a, varargin)
if numel(varargin) > 0
    line_spec = varargin{1};
else
    line_spec = '-';
end
if numel(varargin) > 1
    r_axis_limits = varargin{2};
else
    r_axis_limits = [0, 60];
end
f_ = figure;
a_ = polaraxes;
p_y_ = polarplot([-4, -3:4]*(pi/4), 100.*[relative_class_counts_y(end), relative_class_counts_y]./sum(relative_class_counts_y), [line_spec, 'b'], 'LineWidth', 2);
hold on
p_a_ = polarplot([-4, -3:4]*(pi/4), 100.*[relative_class_counts_a(end), relative_class_counts_a]./sum(relative_class_counts_a), [line_spec, 'r'], 'LineWidth', 2);
a_.RAxis.TickLabelFormat = '%g%%';
a_.RAxis.Limits = r_axis_limits;
f_.Units = 'inches';
f_.Position = [2, 2, 3, 2.5];
a_.ThetaAxis.TickValues = 0:45:315;
a_.ThetaAxis.TickLabelFormat = '%g\x00B0';
a_.ThetaAxis.Label.String = 'Relative angle';
a_.ThetaAxis.Label.Units = 'inches';
a_.ThetaAxis.Label.Position = [2.2, 1.1];
% title('All neurons');
set_plot_3_4_2021();
[h_, p_] = fishertest([...
    relative_class_counts_y(4), sum(relative_class_counts_y) - relative_class_counts_y(4); ...
    relative_class_counts_a(4), sum(relative_class_counts_a) - relative_class_counts_a(4)])
end