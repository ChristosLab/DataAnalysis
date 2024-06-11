function set_plot_6_16(in_cwt, subj_to_plot, kernel_to_plot, task_i, pow_scale)
    cwt_1 = in_cwt{subj_to_plot, 1, kernel_to_plot};
    cwt_2 = in_cwt{subj_to_plot, 2, kernel_to_plot};
    t_lim = [-1.49, 2.49];
    ratio_range = [.2, 5];
    fig_size = [5 1 4.5 3];
    figure('Units', 'inches','Position', fig_size);
    hold on
    t_cwt_1 = [-1.49, -1.01]; % Time for cue onset
    t_cwt_2 = [-0.99, 2.49]; % Time for cue onset
    if pow_scale == 1    
        cwt_y = ratio_range;
        y_label_to_plot = 'Relative power (linear average)';
        to_plot_1 = nanmean(cwt_1, 3);
        to_plot_2 = nanmean(cwt_2, 3);
    elseif pow_scale == 2
        cwt_y = 10*log10(ratio_range);
        y_label_to_plot = 'Relative power (dB)';
        to_plot_1 = nanmean(10*log10(cwt_1), 3);
        to_plot_2 = nanmean(10*log10(cwt_2), 3);
    end

    imagesc(t_cwt_1, [0, 100], to_plot_1, cwt_y);
    imagesc(t_cwt_2, [0, 100], to_plot_2, cwt_y);
    set(gca, 'YDir','normal'); 
    xlim(t_lim);
    ylim([0, 100]);
    
    ps_segment_steps = [-1.5, -1, 0, 0.5, 2];
    ps_xticks = [-1.5, -1.25, -1, -0.5, 0, 0.25, 0.5, 1.25, 2];
    ps_xlabels = {...
        '-1.5', 'baseline', '-1', 'prep', '0', ...
        'cue', '.5', 'delay', '2'...
        };
    as_segment_steps = [-1.5, -1, 0, 0.1];
    as_xticks = [-1.5, -1.25, -1, -0.5, 0];
    as_xlabels = {...
        '-1.5', 'baseline', '-1', 'prep', '0'...
        };
    if task_i == 1
        steps = ps_segment_steps;
        xt = ps_xticks;
        x_lab = ps_xlabels;
    else
        steps = as_segment_steps;
        xt = as_xticks;
        x_lab = as_xlabels;
    end
    for i = 1:numel(steps)
        line([0,0] + steps(i), [1, 101], 'Color', 'k', 'LineStyle', '--');    
    end
    h = colorbar();
    xticks(xt);
    xticklabels(x_lab);
    xlabel('Time from cue onset (s)')
    ylabel('Frequency (Hz)');
    ylabel(h, y_label_to_plot);
    % 
    set(gca, 'FontSize', 9);
    set(gca, 'FontWeight', 'bold');
    end