function neuron_tuning_curve_plot(inputs, best_class_cell, labels, colormap)
%TUNING_CURVE_PLOT plots the tuning curves for 8-class responses
%   inputs - cells where each element is an n-by-8 array
font_size = 10;
fig_size = [2.5, 2];
figure('Unit', 'inches', 'Position', [2, 2, fig_size]);
hold on
x_vals = [0:45:360] - 180;
for i = 1:numel(inputs)
    current_inputs = zeros(size(inputs{i}) + [0, 1]); % Add a duplicate for symmetry
    for j = 1:numel(best_class_cell{i})
        %   So that indices increase anti-clockwise
        current_class_order_ =  get_relative_rad(best_class_cell{i}(j), [8:-1:1, 8])+ 4;
        current_inputs(j, :) = inputs{i}(j, current_class_order_);
    end
    current_mean = mean(current_inputs, 1);
    current_sem = std(current_inputs, [], 1)/sqrt(size(current_inputs, 1));
    l(i) = plot(x_vals, current_mean, 'Color', colormap{i}, 'DisplayName', [labels{i}, sprintf(' (N = %d)', size(current_inputs, 1))]);
    fill([x_vals, flip(x_vals)], [current_mean + current_sem, fliplr(current_mean -current_sem)], colormap{i}, 'FaceAlpha', 0.15, 'EdgeAlpha', 0.15)
end
ylabel('Evoked firing rate (sp/s)')
xlabel('Relative to preferred cue (rad)')
xticks(-180:90:180)
xticklabels({'-\pi', '-\pi/2' ,'0' ,'\pi/2', '\pi'});
ylim([-2, 15]);
legend(l)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',font_size,'FontWeight','Bold', 'LineWidth', 1);
end