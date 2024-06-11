function zw_set_plot_5_5(flag)
xlim([-1,.5])
hold on
plot([0, 0], [0, 100], '--w', 'LineWidth', 2)
xlabel('Time from cue onset (s)');
ylabel('Frequency (Hz)');
if flag == 3
    colormap(hot)
    h = colorbar();
    ylabel(h, 'F')
elseif flag == 1
    colormap()
    h = colorbar();
    ylabel(h, 'Adult - Young');
    h.Ticks = [-1 ,0 ,1];
    h.TickLabels = {'-', 'n.s.', '+'};

elseif flag == 2
    colormap()
    h = colorbar();
    ylabel(h, 'AS - PS');
    h.Ticks = [-1 ,0 ,1];
    h.TickLabels = {'-', 'n.s.', '+'};
elseif flag == 4
    h = colorbar();
    ylabel(h, '\DeltaP_a_d_u_l_t - \DeltaP_y_o_u_n_g')
end
h.FontSize = 18;
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');
end