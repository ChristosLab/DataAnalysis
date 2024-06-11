%%  Path setting
project_dir = 'F:\CCLab\Projects\LFP_project_v0.2';
code_lib = 'Wang_Code';
fig_lib = 'Wang_Figure';
output_database = 'Wang_Database';
antisac_database = fullfile(project_dir, '//LFPDatabychannel_Antisac//'); % Repository to antisaccade data.
prosac_database = fullfile(project_dir, '//LFPDatabychannel_Prosac//'); % Repository to prosaccade data.
file_identifier = '*CH*.mat';
antisac_files = dir(fullfile(antisac_database, file_identifier));
prosac_files = dir(fullfile(prosac_database, file_identifier));
all_files = [antisac_files; prosac_files];
addpath(fullfile(project_dir, code_lib));
%%  Load LFP_repo
load(fullfile(project_dir, output_database, 'LFP_repo_w_pretrial.mat'));
%%  Define target labels and values for entry extraction
channels = unique([LFP_repo.channel_id]);
tasks = unique([LFP_repo.task_id]);
sessions = unique([LFP_repo.session_id]);
monkeys = unique([LFP_repo.monkey_id]);
groups = cell([numel(channels), numel(tasks), numel(sessions), numel(monkeys)]);
entry_labels = {'channel_id', 'task_id', 'session_id', 'monkey_id'};
%%  Extract rows of LFP_repo according to entry_labels
% field_names = fieldnames(LFP_repo);
% LFP_repo_lite = rmfield(LFP_repo, field_names(~ismember(field_names, entry_labels)));

LFP_repo_lite = [];
for e_l = entry_labels
    LFP_repo_lite = [LFP_repo_lite, [LFP_repo.(e_l{1})]'];
end
%%  Group trial ids by the labels and save in cell arrays
tic
for a = 1:size(groups, 1)
    for b = 1:size(groups, 2)
        for c = 1:size(groups, 3)
            c
            for d = 1:size(groups, 4)
                groups{a, b, c, d} = ismember(...
                    LFP_repo_lite, [channels(a), tasks(b), sessions(c), monkeys(d)], 'rows'...
                    );
%                 groups{a, b, c, d} = zw_extract_entry(...
%                     LFP_repo_lite, ...
%                     entry_labels, ...
%                     {channels(a), tasks(b), sessions(c), monkeys(d)}...
%                     );
            end
            toc
        end
    end
end
%%
save(fullfile(project_dir, output_database, 'groups_5_5.mat'), 'groups');
% save(fullfile(project_dir, output_database, 'groups_4_29.mat'), 'groups');
%%  Group periodogram parameters
dur = 1; %  Length of target signal;
lb = 1;
ub = 60;
f_interval = 2;
n_band = floor((ub - lb)/f_interval) + 1;
fs = 500;
%%  Saves average periodogram /channel/task/session/subject in a 5-D matrix
%   5th dimension being time points.
%   THIS NEEDS TO BE A FUNCTION. ZW, May 4th 2020.
periodos = zeros([size(groups), n_band]);
tic
for a = 1:size(groups, 1)
    for b = 1:size(groups, 2)
        for c = 1:size(groups, 3)
            c
            for d = 1:size(groups, 4)
                if ~any(groups{a, b, c ,d})
                    periodos(a, b, c, d, :) = nan;
                    continue;
                end 
                periodos(a, b, c, d, :) = zw_group_power(...
                    periodos(a, b, c, d, :), LFP_repo(groups{a, b, c ,d}), ...
                    dur, lb, ub, f_interval, fs...
                    );
            end
            toc
        end
    end
end
%%
save(fullfile(project_dir, output_database, 'periodos_5_5.mat'), 'periodos');
% save(fullfile(project_dir, output_database, 'periodos_4_29.mat'), 'periodos');
%%
last_young_session = [20, 26, 66, 27];
%   Convert  dummy code
[~, last_young_session] = max(last_young_session == sessions', [], 1);
%%
periodos_young_ps = [];
periodos_young_as = [];
periodos_adult_ps = [];
periodos_adult_as = [];
tic
for a = 1:size(groups, 1)
    for b = 1:size(groups, 2)
        for c = 1:size(groups, 3)
            c
            for d = 1:size(groups, 4)
                p_ = periodos(a, b, c, d, :);
                p_ = reshape(p_, 1, numel(p_));
                if any(isnan(p_))
                    continue
                end
                if c > last_young_session(d)
                    if b == 1
                        periodos_adult_ps = [periodos_adult_ps; p_];
                    else
                        periodos_adult_as = [periodos_adult_as; p_];
                    end
                else
                    if b == 1
                        periodos_young_ps = [periodos_young_ps; p_];
                    else
                        periodos_young_as = [periodos_young_as; p_];
                    end
                end
            end
            toc
        end
    end
end
%%
save(...
    fullfile(project_dir, output_database, 'periodos_by_stage_by_task_5_5.mat'), ...
    'periodos_adult_ps', 'periodos_adult_as', 'periodos_young_ps', 'periodos_young_as'...
    );
% save(...
%     fullfile(project_dir, output_database, 'periodos_by_stage_by_task_4_29.mat'), ...
%     'periodos_adult_ps', 'periodos_adult_as', 'periodos_young_ps', 'periodos_young_as'...
%     );
%%  Plot periodogram by task-by-stage
f_for_plot =( lb:f_interval:ub) + f_interval/2;
fig_size = [5 1 8 5];

%%  yp-ya
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram(periodos_young_ps, f_for_plot, [.0, .0, .9]);
lines(2) = zw_plot_periodogram(periodos_young_as, f_for_plot, [.1, .1, .4]);
legend(lines, {'Young-PS', 'Young-AS'});
print(gcf, fullfile(project_dir, fig_lib, 'yp-ya'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib, 'yp-ya'), '-dtiff', '-r400')

%%  ap-aa
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram(periodos_adult_ps, f_for_plot, [.9, .0, .0]);
lines(2) = zw_plot_periodogram(periodos_adult_as, f_for_plot, [.4, .1, .1]);
legend(lines, {'Adult-PS', 'Adult-AS'});
print(gcf, fullfile(project_dir, fig_lib, 'ap-aa'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib, 'ap-aa'), '-dtiff', '-r400')
%%  all
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram(periodos_adult_ps, f_for_plot, [.9, .0, .0]);
lines(2) = zw_plot_periodogram(periodos_adult_as, f_for_plot, [.4, .1, .1]);
lines(3) = zw_plot_periodogram(periodos_young_ps, f_for_plot, [.0, .0, .9]);
lines(4) = zw_plot_periodogram(periodos_young_as, f_for_plot, [.1, .1, .4]);
significant = f_for_plot(f1(:, 2) < .05);
for i = 1:length(significant)
    lines(5) = plot(significant(i) + [-f_interval/2 ,0, f_interval/2], zeros(1, 3) + 0.01, 'Color', [.6, .6, .6], 'LineWidth', 4);
end
legend(lines, {'Adult-PS', 'Adult-AS', 'Young-PS', 'Young-AS', 'Significant stage main effect'});
xlim([0, 60])
print(gcf, fullfile(project_dir, fig_lib, 'all'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib, 'all'), '-dtiff', '-r400')
%%  ap-yp
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram(periodos_adult_ps, f_for_plot, [.9, .0, .0]);
lines(2) = zw_plot_periodogram(periodos_young_ps, f_for_plot, [.0, .0, .9]);
legend(lines, {'Adult-PS', 'Young-PS'});
print(gcf, fullfile(project_dir, fig_lib, 'ap-yp'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib, 'ap-yp'), '-dtiff', '-r400')
%%  aa-ya
figure('Units', 'inches','Position', fig_size);
hold on
lines = [];
lines(1) = zw_plot_periodogram(periodos_adult_as, f_for_plot, [.4, .1, .1]);
lines(2) = zw_plot_periodogram(periodos_young_as, f_for_plot, [.1, .1, .4]);
legend(lines, {'Adult-AS', 'Young-AS'});
print(gcf, fullfile(project_dir, fig_lib, 'aa-ya'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib, 'aa-ya'), '-dtiff', '-r400')
%%
f1 = zeros(size(f_for_plot, 2), 2);
f2 = zeros(size(f_for_plot, 2), 2);
f12 = zeros(size(f_for_plot, 2), 2);
for i = 1:size(f_for_plot, 2)
    tbl = zw_quick_anova(periodos_young_ps, periodos_young_as, periodos_adult_ps, periodos_adult_as, i);
    f1(i, :) = [tbl{2, 6}, tbl{2, 7}];
    f2(i, :) = [tbl{3, 6}, tbl{3, 7}];
    f12(i, :) = [tbl{4, 6}, tbl{4, 7}];
end
%%
figure('Units', 'inches','Position', fig_size)
plot(f_for_plot, f1(:, 2)', 'b', 'LineWidth', 3)
hold on
plot(f_for_plot, f2(:, 2)', 'r', 'LineWidth', 3)
hold on
plot(f_for_plot, f12(:, 2)', 'g', 'LineWidth', 3)
hold on
plot([f_for_plot(1), f_for_plot(end)], [0.05, 0.05], 'LineStyle', '--', 'Color', [.5,.5,.5], 'LineWidth', 3)
legend({'Stage', 'Task', 'Interaction'})
xlabel('Passband central frequency (Hz)');
ylabel('2-way ANOVA p value');
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');
print(gcf, fullfile(project_dir, fig_lib, 'p_val'), '-dpng', '-r400')
print(gcf, fullfile(project_dir, fig_lib, 'p_val'), '-dtiff', '-r400')