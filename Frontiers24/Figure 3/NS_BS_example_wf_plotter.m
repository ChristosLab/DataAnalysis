clear all; close all; clc;
load('BS_wf.mat')

%%
fig=figure;
time_f=[-11:29]*25;
% p1=plot(time_f,x.ns_wf, '*', 'Color',[0.4 0 0], 'MarkerSize',3); hold on
p2=plot(time_f,x.bs_wf, '*', 'Color',[0 0 0.4], 'MarkerSize',3); hold on

%%
time_c=time_f;
% time_c=[3:29]*25;
spline_points = 100;  % Number of points on the spline curve
spline_time = linspace(min(time_c), max(time_c), spline_points);

%%
% fig1=figure;
% plot(time,x.ns_wf, 'o', 'Color',[0.4 0 0]); hold on
% ns_spline_values = spline(time_c, x.ns_fit, spline_time);
% [max_value, max_idx] = max(ns_spline_values);
% p3=plot(spline_time, ns_spline_values, '-', 'Color', 'r', 'LineWidth', 2);
% plot(spline_time(max_idx),max_value, 'o', 'Color','r', 'MarkerSize',10);
% plot(time,x.ns_fit, '-', 'Color','r','LineWidth',2);

% set(gca, 'FontSize', 15);
% xlabel('Waveform Duration ($\mu$s)', 'Interpreter', 'latex', 'FontSize',15); 
% ylabel('Normalized amplitude', 'Interpreter', 'latex', 'FontSize',15);
% title(sprintf('Narrow Unit'), 'Interpreter', 'latex', 'FontSize',15)
% box off

hold on
%%
% fig2=figure;
% plot(time,x.bs_wf, 'o', 'Color',[0 0 0.4]); hold on
bs_spline_values = spline(time_c, x.bs_fit, spline_time);
[max_value, max_idx] = max(bs_spline_values(12:end));
p4=plot(spline_time, bs_spline_values, '-', 'Color', 'b', 'LineWidth', 2);
plot(spline_time(max_idx+11),max_value, 'o', 'Color',[0 0 0.4], 'MarkerSize',10);
% plot(time,x.bs_fit, '-', 'Color','b','LineWidth',2);

ylim([-1.2 0.8])
set(gca, 'FontSize', 15);
xlabel('Waveform Duration ($\mu$s)', 'Interpreter', 'latex', 'FontSize',15); 
ylabel('Normalized amplitude', 'Interpreter', 'latex', 'FontSize',15);
% title(sprintf('Broad unit'), 'Interpreter', 'latex', 'FontSize',15)

% legend([p1, p3], 'Narrow Empirical', 'Narrow Fitted')
legend([p2, p4], 'Broad Empirical', 'Broad Fitted')
box off

%%
set(fig, 'Renderer', 'painters');
exportgraphics(fig, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure1\BS_full_wf.emf', 'ContentType', 'vector');
print(fig, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure1\BS_full_wf.pdf', '-vector', '-bestfit', '-dwinc');













