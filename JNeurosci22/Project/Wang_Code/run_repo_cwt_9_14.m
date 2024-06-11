clear
%%
zw_setpath
%%
addpath(fullfile(project_dir, code_lib, 'External\fieldtrip\'));
ft_defaults
addpath(genpath(fullfile(project_dir, code_lib, 'External\chronux_2_12\')));
%%
load(fullfile(project_dir, output_database, 'LFP_repo_lite.mat'));
%%  CWT parameters
pre_dur = 0.5;
post_dur = 3.5;
normalizer = 3;
kernel_flag = 3;
down_sample = 10;
fs = 500;
f_range = 2:2:100;
avg_method = 1; %   0 - arithmetic; 1 - geometric
%%  Compute single trial post cwt
seg_size = 10000;
for i = 1:ceil(size(LFP_repo_lite, 1)/seg_size)
    tt = clock;
    seg_start = 1 + (i - 1)*seg_size;
    seg_end = min(seg_start + seg_size -1, size(LFP_repo_lite, 1));
    fprintf(...
        'Segment %d to %d started at %s\n', ...
        seg_start, seg_end, string(datetime)...
        );
    repo_cwt = zw_repo_cwt_no_pre(...
        LFP_repo(1, seg_start:seg_end), ...
        post_dur, normalizer, ...
        kernel_flag, down_sample, f_range, fs...
        );
    save(...
        fullfile(...
        project_dir, ...
        output_database, ...
        sprintf('repo_cwt_post_only_%d_to_%d.mat', seg_start, seg_end)...
        ), ...
        'repo_cwt'...
        );
    repo_cwt = [];
end
%%  Load and recreate complete cwt structure
repo_cwt_total = struct('cwt', []);
seg_size = 10000;
for i = 1:ceil(size(LFP_repo_lite, 1)/seg_size)
    seg_start = 1 + (i - 1)*seg_size;
    seg_end = min(seg_start + seg_size -1, size(LFP_repo_lite, 1));
    fprintf(...
        'Loading saved cwt from segment %d to %d started at %s\n', ...
        seg_start, seg_end, string(datetime)...
        );
    load(...
        fullfile(...
        project_dir, ...
        output_database, ...
        sprintf('repo_cwt_post_only_%d_to_%d.mat', seg_start, seg_end)...
        )...
        );
    repo_cwt_total = [repo_cwt_total, repo_cwt];
    repo_cwt = [];
end
repo_cwt_total = repo_cwt_total(2:end);
%%  Compute TFR patch averages
fr_fs = 0.5;
epoch_fs = 50;
target_frs = [1, 2; 3, 4; 5, 6; 5, 9; 7, 15; 16, 50];
target_epochs = [1, 50; 51, 75; 76, 150; 51, 55];
target_output = zeros([...
    size(target_frs, 1), ...
    size(target_epochs, 1), ...
    numel(repo_cwt_total)...
    ]);
for i = 1:size(target_frs, 1)
    for j = 1:size(target_epochs, 1)
        for k = 1:numel(repo_cwt_total)
            target_output(i, j, k) = zw_tfr_patch(...
                target_frs(i, :), ...
                target_epochs(j, :), ...
                repo_cwt_total(k).cwt...
                );
        end
    end
end
clear repo_cwt_total
%%  Save target output
%   Output is a band-by-epoch-by-trial tensor.
dt = datetime;
dt.Format = 'u_MM_dd_HH_mm_ss';
save(fullfile(project_dir, output_database, sprintf('target_output_%s.mat', string(dt))), 'target_output');
%%  1-way ANOVA examining spatial tuning for each recording site
analysis_flag = 0;
%   Load group information
load(fullfile(project_dir, output_database, 'groups_9_14.mat'));
%%
%  Load target summary output
file_ = uigetfile(fullfile(project_dir, output_database, 'target_output*'), 'Please select target output file: ');
load(fullfile(project_dir, output_database, file_));
%%
%   Step1: Do the ANOVA on each entry (Subject X task X session X channel);
%   Assign to each entry a label, store label info in a 2-D array
analysis_output = struct('tb', [], 'stats', []);
analysis_output_groups = zeros(0, 4);
combine_task_variate_flag = 1;
site_counter = 0;
tic
for d = 1:size(groups, 4) % 'monkey_id'
    for b = 1:size(groups, 2) % 'task_id'
        for c = 1:size(groups, 3) % 'session_id'
            for a = 1:size(groups, 1) % 'channel_id'
                [referenced_cells, referenced_subs] = zw_cell_reference([a, b, c, d, 0], groups);
                if combine_task_variate_flag
                    [referenced_cells, referenced_subs] = combine_task_variate(referenced_cells, referenced_subs);
                end
                if numel(referenced_cells) ~= 0
                    site_counter = site_counter + 1;
                    for fr = 1:size(target_output, 1)
                        for ep = 1:size(target_output, 2)
                            target_output_single = squeeze(target_output(fr,ep,:));
                            [temp_tb, temp_stats] = analysis_call(referenced_cells, referenced_subs, target_output_single, analysis_flag);
                            analysis_output(fr, ep, site_counter).tb = temp_tb;
                            analysis_output(fr, ep, site_counter).stats = temp_stats;
                        end
                    end
                    analysis_output_groups = [analysis_output_groups; [d, b, c ,a]]; 
                    toc
                end
            end
        end
    end
end
toc
%%  Save analysis_output
dt = datetime;
dt.Format = 'u_MM_dd_HH_mm_ss';
save(fullfile(project_dir, output_database, sprintf('analysis_output_%s.mat', string(dt))), 'analysis_output', 'analysis_output_groups');
%%  Load analysis_output
file_ = uigetfile(fullfile(project_dir, output_database, 'analysis_output*'), 'Please select analysis output file: ');
load(fullfile(project_dir, output_database, file_));
%%  Add an aged/young binary column called "aged"
%   0 - young, 1 - adult
last_young_session = [20, 26, 66, 27];
aged = zeros(size(analysis_output_groups, 1), numel(last_young_session));
for i = 1:numel(last_young_session)
    aged(:, i) = [analysis_output_groups(:, 1) == i].*[analysis_output_groups(:, 3) > last_young_session(i)];
end
aged = sum(aged, 2);
%%  Y-P
%                 subject, task, session, channel, stage
sub_select_rule = [0,    1, 0,    0,    0; ...
    1000, 1, 1000, 1000, 0]; %  First row: Lower bound (closed)
%  Second row: Upper
selected_ = and(prod([analysis_output_groups, aged] >= sub_select_rule(1, :), 2)...
    , prod([analysis_output_groups, aged] <= sub_select_rule(2, :), 2));
%%  A-P
sub_select_rule = [0,    1, 0,    0,    1; ...
    1000, 1, 1000, 1000, 1]; %  First row: Lower bound (closed)
%                                                       Second row: Upper
selected_ = and(prod([analysis_output_groups, aged] >= sub_select_rule(1, :), 2)...
    , prod([analysis_output_groups, aged] <= sub_select_rule(2, :), 2));
%%  Y-A
%                 subject, task, session, channel, stage
sub_select_rule = [0,    2, 0,    0,    0; ...
    1000, 2, 1000, 1000, 0]; %  First row: Lower bound (closed)
%                                                       Second row: Upper
selected_ = and(prod([analysis_output_groups, aged] >= sub_select_rule(1, :), 2)...
    , prod([analysis_output_groups, aged] <= sub_select_rule(2, :), 2));
%%  A-A
%                 subject, task, session, channel, stage
sub_select_rule = [0,    2, 0,    0,    1; ...
    1000, 2, 1000, 1000, 1]; %  First row: Lower bound (closed)
%                                                       Second row: Upper
selected_ = and(prod([analysis_output_groups, aged] >= sub_select_rule(1, :), 2)...
    , prod([analysis_output_groups, aged] <= sub_select_rule(2, :), 2));
%%  Plot 1-way ANOVA p values in a histogram
band_names = {'Delta', 'Theta', 'Alpha', 'Alpha-Beta', 'Beta', 'Gamma'};
epoch_names = {'Fixation', 'Cue', 'Delay', 'Cue'};
stage_names = {'Adolescent', 'Adult'};
task_names = {'Pro-Saccade', 'Anti-Saccade'};
face_color = [0.05, 0.4, 0.95; 0.8, 0, 0.1];

%
task_idx = sub_select_rule(1, 2);
stage_idx = sub_select_rule(1, end) + 1;
stage_counts = [sum(aged == 0), sum(aged == 1)];
stage_x_task_counts = zeros(2, 2);
for i = 1:2
    for j = 1:2
        stage_x_task_counts(i, j) = sum((aged == (i - 1)).*analysis_output_groups(:, 2) == j);
    end
end
%
yl = [0 0.4];
alpha_level = 0.05;
if task_idx == 1
    j_range = 1:3;
else
    j_range = [1,4];
end
%
figure('Unit', 'inches', 'Position', [0, 0, 2*size(target_frs, 1), 1 + 2*numel(j_range)])
for i = 1:size(target_frs, 1)
    for j= 1:numel(j_range)
        subplot(numel(j_range), size(target_frs, 1), i + size(target_frs, 1)*(j - 1))
        hold on
        histogram(...
            ref_nested_cell([2, 6], ...
            {analysis_output(i, j_range(j), selected_).tb}), ...
            'BinWidth', alpha_level, ...
            'Normalization', 'probability', ...
            'FaceColor', face_color(task_idx, :));
        line([0, 0] + alpha_level, yl, 'LineStyle', '--', 'Color', [0, 0.6, 0])
        ylim(yl)
        if and(i == 1, j == ceil(numel(j_range)/2))
            ylabel(sprintf('%s\nProportion LFP units', epoch_names{j_range(j)}))
        elseif i == 1
            ylabel(sprintf('%s\n ', epoch_names{j_range(j)}))
        elseif and(j == numel(j_range), i == ceil(size(target_frs, 1)/2))
            xlabel('Stimulus class: 1-way ANOVA p-value')
        end
        if j == 1
            title(sprintf('%s(%d-%d Hz)', band_names{i}, target_frs(i, 1)/fr_fs, target_frs(i, 2)/fr_fs))
        end
        if and(i == 1, j == 1)
            %         text(alpha_level, yl(2)*0.8, sprintf('p value = %.2f', alpha_level))
            text(alpha_level, yl(2)*0.8, num2str(alpha_level))
        end
    end
end
fig_name = sprintf('Class Tuning of %s LFP units in %s (N = %d)', stage_names{stage_idx}, task_names{task_idx}, stage_x_task_counts(stage_idx, task_idx));
sgtitle(fig_name,'FontWeight', 'Bold', 'FontName','Times New Roman','FontSize',18);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Refer back to LFP sites based on analysis stats
%  Y-P
%                 subject, task, session, channel, stage
sub_select_rule = [0,    1, 0,    0,    0; ...
    1000, 1, 1000, 1000, 0]; %  First row: Lower bound (closed)
%  Second row: Upper
selected_yp = and(prod([analysis_output_groups, aged] >= sub_select_rule(1, :), 2)...
    , prod([analysis_output_groups, aged] <= sub_select_rule(2, :), 2));
%  A-P
sub_select_rule = [0,    1, 0,    0,    1; ...
    1000, 1, 1000, 1000, 1]; %  First row: Lower bound (closed)
%                                                       Second row: Upper
selected_ap = and(prod([analysis_output_groups, aged] >= sub_select_rule(1, :), 2)...
    , prod([analysis_output_groups, aged] <= sub_select_rule(2, :), 2));
%
alpha_level = 0.05;
band_names = {'Delta', 'Theta', 'Alpha', 'Alpha-Beta', 'Beta', 'Gamma'};
epoch_names = {'Fixation', 'Cue', 'Delay', 'Cue'};
stage_names = {'Adolescent', 'Adult'};
task_names = {'Pro-Saccade', 'Anti-Saccade'};
%
task_idx = sub_select_rule(1, 2);
stage_idx = sub_select_rule(1, end) + 1;
stage_counts = [sum(aged == 0), sum(aged == 1)];
% sig_sites = zeros(size(selected_, 1), size(target_frs, 1), numel(j_range));
sig_yp_p_value = cell(size(target_frs, 1), size(target_epochs, 1));
sig_ap_p_value = cell(size(target_frs, 1), size(target_epochs, 1));
for i = 1:size(target_frs, 1)
    for j= 1:size(target_epochs, 1)
        sig_sites_yp(:, i, j) = ismember(1:size(selected_yp, 1), ref_func_output(...
            find(selected_yp), ...
            find(ref_nested_cell([2, 6], {analysis_output(i, j, selected_yp).tb}) < alpha_level)...
            ));
        sig_yp_p_value{i, j} = ref_nested_cell([2, 6], {analysis_output(i, j, sig_sites_yp(:, i, j)).tb});
        sig_sites_ap(:, i, j) = ismember(1:size(selected_ap, 1), ref_func_output(...
            find(selected_ap), ...
            find(ref_nested_cell([2, 6], {analysis_output(i, j, selected_ap).tb}) < alpha_level)...
            ));
        sig_ap_p_value{i, j} = ref_nested_cell([2, 6], {analysis_output(i, j, sig_sites_ap(:, i, j)).tb});
    end
end
%%
crit_fr = 4;
crit_epoch = 3;
%   sig_sites_yp(sites, frs, epochs)
t_yp = find(sig_sites_yp(:,crit_fr,crit_epoch));
t_ap = find(sig_sites_ap(:,crit_fr,crit_epoch));
t_yp_n = find((analysis_output_groups(:, 2) == 1).*(aged == 0).*~sig_sites_yp(:,crit_fr,crit_epoch));
t_ap_n = find((analysis_output_groups(:, 2) == 1).*(aged == 1).*~sig_sites_ap(:,crit_fr,crit_epoch));
%   Anti-saccade trials session-matched to the pro-saccade trials
t_ya = find((analysis_output_groups(:, 2) == 2).*ismember(analysis_output_groups(:, [1,3,4]), analysis_output_groups(t_yp, [1,3,4]), 'rows'));
t_aa = find((analysis_output_groups(:, 2) == 2).*ismember(analysis_output_groups(:, [1,3,4]), analysis_output_groups(t_ap, [1,3,4]), 'rows'));
t_ya_n = find((analysis_output_groups(:, 2) == 2).*ismember(analysis_output_groups(:, [1,3,4]), analysis_output_groups(t_yp_n, [1,3,4]), 'rows'));
t_aa_n = find((analysis_output_groups(:, 2) == 2).*ismember(analysis_output_groups(:, [1,3,4]), analysis_output_groups(t_ap_n, [1,3,4]), 'rows'));
%%  Reference based on t_yp/t_yp_n/t_ap/t_ap_n
stage_idx = 1;
sig_seq_for_plot = max_n_ind(sig_yp_p_value{crit_fr, crit_epoch}, -5);
t_for_plot = t_yp(sig_seq_for_plot)';
% t_for_plot = t_yp_n(70);
%%
stage_idx = 2;
sig_seq_for_plot = max_n_ind(sig_ap_p_value{crit_fr, crit_epoch}, -10);
t_for_plot = t_yp(sig_seq_for_plot)';

%%  Plot tuning of single LFP units
for p = t_for_plot
    % analysis_output_group dimension labels: 'Monkey', 'Task', 'Session', 'channel'
    %                          groups labels: 'channel', 'Task', 'Session', 'Monkey'
    [referenced_cells, referenced_subs] = zw_cell_reference([analysis_output_groups(p, [4,2,3,1]), 0], groups);
    %  Plot for single sessisson all classes
    figure('Unit', 'inches', 'Position', [0, 0, 7, 7]);
    subplot_positions = [6, 3, 2, 1, 4, 7, 8, 9];
    mean_target_output = zeros([size(target_output, 1), ...
        size(target_output, 2), ...
        numel(referenced_subs)]);
    cwt_range = [10e9,10e-10];
%     for i = 1:numel(referenced_subs) %  Classes
%         %   1: Raw TFR plotting
%         mean_cwt = [];
%         cwts_ = {repo_cwt_total(referenced_cells{i}).cwt};
%         patch_tfr = target_output(:, :, referenced_cells{i});
%         for j = 1:numel(cwts_) %    Trials
%             if and(numel(cwts_{j})>1, ~any(isnan(cwts_{j})))
%                 mean_cwt = cat(3, mean_cwt, cwts_{j}); %   Mean of one class acriss trials
%             end
%         end
%         cwt_range(1) = min(cwt_range(1), min(min(mean(10*log10(mean_cwt), 3))));
%         cwt_range(2) = max(cwt_range(2), max(max(mean(10*log10(mean_cwt), 3))));
%     end
    for i = 1:numel(referenced_subs) %  Classes
        %   1: Raw TFR plotting
        mean_cwt = [];
        cwts_ = {repo_cwt_total(referenced_cells{i}).cwt};
        patch_tfr = target_output(:, :, referenced_cells{i});
        for j = 1:numel(cwts_) %    Trials
            if and(numel(cwts_{j})>1, ~any(isnan(cwts_{j})))
                mean_cwt = cat(3, mean_cwt, cwts_{j}); %   Mean of one class acriss trials
            end
        end
        mean_cwt = mean(10*log10(mean_cwt), 3);
        subplot(3, 3, subplot_positions(i));
        imagesc([1/epoch_fs, 175/epoch_fs], [2, 100], mean_cwt, [10*log10(0.0005), 10*log10(0.02)])
%         imagesc([1/epoch_fs, 175/epoch_fs], [2, 100], mean_cwt, cwt_range/2)
        xlim([1.5, 3])
        set(gca, 'YDir', 'normal')
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
        %   2: Summary stats computation
        mean_target_output(:, :, i) = nanmean(patch_tfr, 3);
    end
    % Middle polar-plot of power
    crit_tfr = mean_target_output(crit_fr, crit_epoch, :);
    crit_tfr = reshape(crit_tfr, 1, numel(crit_tfr));
    crit_tfr = crit_tfr/max(crit_tfr);
    crit_tfr = [crit_tfr, crit_tfr(1)];
    subplot(3, 3, 5);
    polarplot(([referenced_subs, referenced_subs(1)] - 1)*(pi/4), crit_tfr);
    title('Proportion power of preferred location')
%     
    fig_name = sprintf('%s LFP unit #%d with significant delay %s tuning', stage_names{stage_idx}, p, band_names{crit_fr});
    sgtitle(fig_name,'FontWeight', 'Bold', 'FontName','Times New Roman','FontSize',14);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%%  Save significant LFP sites and channel/session information into CSV file
sig_sites_p = sig_sites_yp + sig_sites_ap;
lfp_sites_tbl = cell(size(analysis_output_groups, 1), 0);
monkey_id = {'ind', 'jac', 'lem', 'ken'};
task_names = {'ODR', 'Anti-saccade'};
for i =1:size(analysis_output_groups, 1)
    lfp_sites_tbl{i, 1} = monkey_id{analysis_output_groups(i, 1)};
    lfp_sites_tbl{i, 2} = entry_label_values{3}(analysis_output_groups(i, 3));
    lfp_sites_tbl{i, 3} = task_names{analysis_output_groups(i, 2)};
    lfp_sites_tbl{i, 4} = entry_label_values{1}(analysis_output_groups(i, 4));
    lfp_sites_tbl{i, 5} = sig_sites_p(i, 6, 3);
    lfp_sites_tbl{i, 6} = sig_sites_p(i, 4, 3);
end
lfp_sites_tbl = cat(2, num2cell([1:size(analysis_output_groups, 1)]'), lfp_sites_tbl);
lfp_sites_tbl = cell2table(...
    lfp_sites_tbl, ...
    'VariableNames', cat(2, 'No', entry_labels([4, 3, 2, 1]), {'ODR_delay_gamma_tuning', 'ODR_delay_alpha_beta_tuning'})...
    );
writetable(lfp_sites_tbl, fullfile(...
        project_dir, ...
        output_database, ...
        'LFP_sites.csv'...
        ));
%%
writetable(lfp_sites_tbl, fullfile(...
        project_dir, ...
        output_database, ...
        'LFP_sites.csv'...
        ));
%%  significant LFP sites and channel/session information LITE
sig_sites_p = sig_sites_yp + sig_sites_ap;
lfp_sites_tbl = cell(size(analysis_output_groups, 1), 0);
monkey_id = {'ind', 'jac', 'lem', 'ken'};
task_names = {'ODR', 'Anti-saccade'};
for i =1:size(analysis_output_groups, 1)
    lfp_sites_tbl{i, 1} = analysis_output_groups(i, 1); %   Monkey
    lfp_sites_tbl{i, 2} = entry_label_values{3}(analysis_output_groups(i, 3)); %    Session
    lfp_sites_tbl{i, 3} = analysis_output_groups(i, 2); %   Task
    lfp_sites_tbl{i, 4} = entry_label_values{1}(analysis_output_groups(i, 4)); % Channel (do not enumerate)
    lfp_sites_tbl{i, 5} = sig_sites_p(i, 6, 3);
    lfp_sites_tbl{i, 6} = sig_sites_p(i, 5, 3);
    lfp_sites_tbl{i, 7} = sig_sites_p(i, 4, 3);
    lfp_sites_tbl{i, 8} = sig_sites_p(i, 6, 2);
    lfp_sites_tbl{i, 9} = sig_sites_p(i, 5, 2);
    lfp_sites_tbl{i, 10} = sig_sites_p(i, 4, 2);

end
lfp_sites_tbl = cat(2, num2cell([1:size(analysis_output_groups, 1)]'), lfp_sites_tbl);
lfp_sites_tbl = cell2table(...
    lfp_sites_tbl, ...
    'VariableNames', cat(2, 'No', entry_labels([4, 3, 2, 1]), {'delay_gamma', 'delay_beta', 'delay_alphabeta', 'cue_gamma', 'cue_beta', 'cue_alphabeta'})...
    );
%%  Save LFP sites table as MAT file
fname_ = 'lfp_sites_tbl.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_sites_tbl');
%%  Load LFP sites table as MAT file
fname_ = 'lfp_sites_tbl.mat';
load(fullfile(project_dir, output_database, fname_));
%%  Load neuron table
fname_ = 'neuron_tbl.mat';
load(fullfile(project_dir, output_database, fname_));
%%  LFP - Neuron pairing 
% Needs to switch channel_id from 1-based indexing to actual channel number
neuron_info_columns = 9:14;
% Last column records the number of isolated neurons at site
neuron_info_by_site_ = zeros(size(analysis_output_groups, 1), numel(neuron_info_columns) + 1);
for i = 1:size(analysis_output_groups, 1)
    monkey_id_ = analysis_output_groups(i, 1);
    task_id_ = analysis_output_groups(i, 2);
    session_id_ = analysis_output_groups(i, 3);
    channel_id_ = entry_label_values{1}(analysis_output_groups(i, 4));
    matches_ = ismember(table2array(neuron_tbl(:, [3, 6, 7, 8])), [channel_id_, monkey_id_, session_id_, task_id_], 'row');
    neuron_info_by_site_(i, :) = [any(table2array(neuron_tbl(matches_, neuron_info_columns)), 1), sum(matches_)];
end
for i = 1:numel(neuron_info_columns)
    var_to_add_ = neuron_info_by_site_(:, i);
    lfp_sites_tbl = addvars(lfp_sites_tbl, var_to_add_, 'NewVariableNames', neuron_tbl.Properties.VariableNames(neuron_info_columns(i)));
end
var_to_add_ = neuron_info_by_site_(:, end);
lfp_sites_tbl = addvars(lfp_sites_tbl, var_to_add_, 'NewVariableNames', {'n_neuron'});
var_to_add_ = any(neuron_info_by_site_(:, 1:4) , 2); %    Neuronal response to any of the epochs
lfp_sites_tbl = addvars(lfp_sites_tbl, var_to_add_, 'NewVariableNames', {'responsive'});
lfp_sites_tbl = addvars(lfp_sites_tbl, aged, 'NewVariableNames', {'aged'});
%%  Save neuron-LFP sites table as MAT file
fname_ = 'neuron-lfp_sites_tbl.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_sites_tbl');
%%  Filter out needed neurons
needed_neuron = zeros(size(neuron_tbl, 1), 1);
to_match_ = [entry_label_values{1}(analysis_output_groups(:, 4))', analysis_output_groups(:, [1, 3, 2])];
for i = 1:size(neuron_tbl, 1)
    if any(ismember(to_match_, table2array(neuron_tbl(i, [3, 6, 7, 8])) , 'row'))
        needed_neuron(i) = 1;
    end
end
%%  Save neuron data in one repo
% file_names_ = strcat(neuron_tbl.Filename(find(needed_neuron)), '_', num2str(neuron_tbl.Neuron(find(needed_neuron))));
file_names_ = strcat(neuron_tbl.Filename, '_', num2str(neuron_tbl.Neuron));
for i = 1:numel(file_names_)
    if needed_neuron(i)
        try
            load(fullfile(neuron_database, file_names_{i}))
            neuron_repo(i) = MatData;
            clear MatData
        catch
            needed_neuron(i) = 0;
        end
        i
    end
end
fname_ = 'neuron_repo';
save(fullfile(project_dir, output_database, fname_), 'neuron_repo');
%%  Check for non-exising neuron files
m = 0;
for i =1:numel(neuron_repo)
    if numel(neuron_repo(i).class) < 1
        m = m + 1
    end
end
%%  compute psth either cue or saccade aligned
bin_width = 0.05;  % 50 milliseconds bin
bin_edges_cue = -1:bin_width:3; %   [-1, 3] around cue onset, should work for both tasks
bin_edges_sac = -3:bin_width:1; %   [-5, 2] around sccade onset
t_bins_cue = bin_edges_cue(1:(end-1)) + bin_width/2;
t_bins_sac = bin_edges_sac(1:(end-1)) + bin_width/2;
n_bins_cue = numel(histcounts([], bin_edges_cue));
n_bins_sac = numel(histcounts([], bin_edges_sac));
for i = find(needed_neuron)'
    for j = 1:numel(neuron_repo(i).class)
        psth_cue_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_cue);
        psth_sac_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_sac);
        for m = 1:numel(neuron_repo(i).class(j).ntr)
            TS_cue_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Cue_onT;
            psth_cue_(m, :) = histcounts(TS_cue_, bin_edges_cue)/bin_width;
            if ~isempty(neuron_repo(i).class(j).ntr(m).Saccade_onT)
                TS_sac_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Saccade_onT;
                psth_sac_(m, :) = histcounts(TS_sac_, bin_edges_sac)/bin_width;
            else
                psth_sac_(m, :) = NaN;
            end
        end
        neuron_repo(i).class(j).psth_cue = psth_cue_;
        neuron_repo(i).class(j).psth_sac = psth_sac_;
    end
end
%%
fname_ = 'neuron_repo_w_psth';
save(fullfile(project_dir, output_database, fname_), 'neuron_repo');
%%  neuron class ANOVA
% neuron_analysis_output_cue = struct([]);
neuron_analysis_output_sac = struct([]);
for i = 1:numel(neuron_repo)
    if ~needed_neuron(i) %  Omit neuron if not included in LFP sites
    else
        tic
        n_class_ = numel(neuron_repo(i).class);
%         psth_cue_ = {neuron_repo(i).class.psth_cue};
        psth_sac_ = {neuron_repo(i).class.psth_sac};
        for j = 1:n_bins_cue
%             anova_y_cue_ = [];
%             anova_label_cue_ = [];
            anova_y_sac_ = [];
            anova_label_sac_ = [];
            for m = 1:n_class_
                if m > 8
                    class_id_ = mod(m, 8);
                else
                    class_id_ = m;
                end
%                 anova_y_cue_ = [anova_y_, psth_cue_{m}(:, j)'];
%                 anova_label_cue_ = [anova_label_, class_id_ + zeros(size(psth_cue_{m}(:, j)'))];
                anova_y_sac_ = [anova_y_, psth_sac_{m}(:, j)'];
                anova_label_sac_ = [anova_label_, class_id_ + zeros(size(psth_sac_{m}(:, j)'))];
            end
%             [~, out_tb_, out_stats_] = anova1(anova_y_, anova_label_, 'off');
%             neuron_analysis_output_cue(i, j).tb = out_tb_;
%             neuron_analysis_output_cue(i, j).stats = out_stats_;
            [~, out_tb_, out_stats_] = anova1(anova_y_sac_, anova_label_sac_, 'off');
            neuron_analysis_output_sac(i, j).tb = out_tb_;
            neuron_analysis_output_sac(i, j).stats = out_stats_;
        end
        toc
    end
end
%%  Save needed neuron indices
fname_ = 'needed_neuron';
save(fullfile(project_dir, output_database, fname_), 'needed_neuron');
%%  Save neuron_analysis_output_cue
fname_ = 'neuron_analysis_output_cue';
save(fullfile(project_dir, output_database, fname_), 'neuron_analysis_output_cue');
%%  Save neuron_analysis_output_sac
fname_ = 'neuron_analysis_output_sac';
save(fullfile(project_dir, output_database, fname_), 'neuron_analysis_output_sac');
%%  Compute PEV
neuron_pev_cue = zeros(size(neuron_analysis_output_cue));
neuron_pev_sac = zeros(size(neuron_analysis_output_sac));
for i = 1:size(neuron_analysis_output_cue, 1)
    i
    if needed_neuron(i)
        for j = 1:size(neuron_analysis_output_cue, 2)
            neuron_pev_cue(i, j) = zw_pev(neuron_analysis_output_cue(i, j).tb);
        end
        for j = 1:size(neuron_analysis_output_sac, 2)
            neuron_pev_sac(i, j) = zw_pev(neuron_analysis_output_sac(i, j).tb);
        end
    else
        neuron_pev_cue(i, :) = NaN;        
        neuron_pev_sac(i, :) = NaN;
    end
end
%%  Save PEV
fname_ = 'neuron_pev_cue';
save(fullfile(project_dir, output_database, fname_), 'neuron_pev_cue');
fname_ = 'neuron_pev_sac';
save(fullfile(project_dir, output_database, fname_), 'neuron_pev_sac');
%%  Added stage and tuning info to neuron_tbl
last_young_session = [20, 26, 66, 27];
to_match_ = [entry_label_values{1}(analysis_output_groups(:, 4))', analysis_output_groups(:, [1, 3, 2])];
% to_match_ = [entry_label_values{1}(analysis_output_groups(:, 4))', analysis_output_groups(:, [1, 3]), ones(size(analysis_output_groups, 1), 1)];
for i = 1:size(neuron_tbl, 1)
    neuron_tbl.aged(i) = neuron_tbl.session_id(i) > last_young_session(neuron_tbl.monkey_id(i));
    if needed_neuron(i)
        %   Match both PS and AS to PS tuning
        site_idx_ = find(ismember(to_match_, [table2array(neuron_tbl(i, [3, 6, 7])), 1] , 'row'));
        neuron_tbl.delay_gamma(i) = lfp_sites_tbl.delay_gamma(site_idx_);
        neuron_tbl.delay_beta(i) = lfp_sites_tbl.delay_beta(site_idx_);
        neuron_tbl.delay_alphabeta(i) = lfp_sites_tbl.delay_alphabeta(site_idx_);        
        neuron_tbl.cue_gamma(i) = lfp_sites_tbl.cue_gamma(site_idx_);
        neuron_tbl.cue_beta(i) = lfp_sites_tbl.cue_beta(site_idx_);
        neuron_tbl.cue_alphabeta(i) = lfp_sites_tbl.cue_alphabeta(site_idx_);        
    else
        neuron_tbl.delay_gamma(i) = NaN;
        neuron_tbl.delay_beta(i) = NaN;
        neuron_tbl.delay_alphabeta(i) = NaN;        
        neuron_tbl.cue_gamma(i) = NaN;
        neuron_tbl.cue_beta(i) = NaN;
        neuron_tbl.cue_alphabeta(i) = NaN;    
    end
end
%%  Plot PEV based on various grouping criteria
% ODR_cue_aligned_neuron_PEV_plot_11_16
% ODR_sac_aligned_neuron_PEV_plot_11_16
%%  PEV for TFR
clear neuron_repo neuron_analysis_output_cue neuron_analysis_output_sac neuron_pev_sac neuron_pev_cue
%%  Compute moment-by-moment mean TFR
fr_fs = 0.5;
epoch_fs = 50;
target_frs = [1, 2; 3, 4; 5, 6; 5, 9; 7, 15; 16, 50];
mbm_epochs = [1:175; 1:175]';
lfp_mbm_output = zeros([...
    size(target_frs, 1), ...
    size(mbm_epochs, 1), ...
    numel(repo_cwt_total)...
    ]);
for i = 1:size(target_frs, 1)
        for k = 1:numel(repo_cwt_total)
            k
            lfp_mbm_output(i, :, k) = zw_band_power(...
                target_frs(i, :), ...
                repo_cwt_total(k).cwt...
                );
        end
end
% clear repo_cwt_total
%%
%   Output is a band-by-moment-by-trial tensor.
dt = datetime;
dt.Format = 'u_MM_dd_HH_mm_ss';
save(fullfile(project_dir, output_database, sprintf('lfp_mbm_output_%s.mat', string(dt))), 'lfp_mbm_output');
%%
%   Step1: Do the ANOVA on each entry (Subject X task X session X channel);
%   Assign to each entry a label, store label info in a 2-D array
lfp_mbm_analysis_output = struct('tb', [], 'stats', []);
lfp_mbm_analysis_output_groups = zeros(0, 4);
combine_task_variate_flag = 1;
site_counter = 0;
tic
% for d = 1:size(groups, 4) % 'monkey_id'
%     for b = 1:size(groups, 2) % 'task_id'
%         for c = 1:size(groups, 3) % 'session_id'
%             for a = 1:size(groups, 1) % 'channel_id'
for i = 1:size(analysis_output_groups, 1)
    d = analysis_output_groups(i, 1);
    b = analysis_output_groups(i, 2);
    c = analysis_output_groups(i, 3);
    a = analysis_output_groups(i, 4);
                [referenced_cells, referenced_subs] = zw_cell_reference([a, b, c, d, 0], groups);
                if combine_task_variate_flag
                    [referenced_cells, referenced_subs] = combine_task_variate(referenced_cells, referenced_subs);
                end
%                 if numel(referenced_cells) ~= 0
%                     site_counter = site_counter + 1;
                    for fr = 1:size(target_output, 1)
                        for ep = 1:size(lfp_mbm_output, 2)
                            target_output_single = squeeze(lfp_mbm_output(fr,ep,:));
                            [temp_tb, temp_stats] = analysis_call(referenced_cells, referenced_subs, target_output_single, analysis_flag);
                            lfp_mbm_analysis_output(fr, ep, i).tb = temp_tb;
                            lfp_mbm_analysis_output(fr, ep, i).stats = temp_stats;
                        end
                    end
%                    lfp_mbm_analysis_output_groups = [lfp_mbm_analysis_output_groups; [d, b, c ,a]]; 
                    toc
%                 end
%             end
%         end
%     end
end
toc
%%  Save lfp_mbm_analysis_output
dt = datetime;
dt.Format = 'u_MM_dd_HH_mm_ss';
save(fullfile(project_dir, output_database, sprintf('lfp_mbm_analysis_output_%s.mat', string(dt))), 'lfp_mbm_analysis_output', 'lfp_mbm_analysis_output_groups');
%%  Compute LFP PEV
lfp_pev_cue = zeros(size(lfp_mbm_analysis_output));
for i = 1:size(lfp_mbm_analysis_output, 3)
    i
    for j = 1:size(lfp_mbm_analysis_output, 1)
        for k = 1:size(lfp_mbm_analysis_output, 2)
            lfp_pev_cue(j, k, i) = zw_pev(lfp_mbm_analysis_output(j, k, i).tb);
        end
    end
end
lfp_pev_cue = permute(lfp_pev_cue, [3, 2, 1]); %  Changed to trial-by-time-by-fr
%%
fname_ = 'lfp_pev_cue';
save(fullfile(project_dir, output_database, fname_), 'lfp_pev_cue');
%%
%%  Find missing session files
to_match_ = [entry_label_values{1}(analysis_output_groups(:, 4))', analysis_output_groups(:, [1, 3, 2])];
for i = 1:size(neuron_tbl, 1)
    if ~any(ismember(to_match_, table2array(neuron_tbl(i, [3, 6, 7, 8])) , 'row'))
        neuron_from_missing_file(i) = 1;
    end
end
%%
missing_files = table(unique(neuron_tbl.Filename(find(neuron_from_missing_file))), 'VariableNames', {'Filename'});
writetable(missing_files, fullfile(...
        project_dir, ...
        output_database, ...
        'Missing_files.csv'...
        ));

%%
function [out_tb, out_stats] = analysis_call(referenced_cells, referenced_subs, target_output, flag)
%   Squeeze inputs to column vectors
switch flag
    case 0 %    1-way ANOVA
        anova_y = [];
        anova_label = [];
        for i = 1:numel(referenced_cells)
            temp_ = target_output([referenced_cells{i}]);
            anova_y = [anova_y; temp_]; %  Make a point to only operate on 1-D arrays; Keep it simple.
            anova_label = [anova_label; ones(numel(temp_), 1)*referenced_subs(i)];
        end
        [~, out_tb, out_stats] = anova1(anova_y, anova_label, 'off');
    otherwise
        disp('Unrecognized analysis type flag\n')
end
end


