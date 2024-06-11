bin_width = 0.05;  % 50 milliseconds bin
bin_edges_cue = -1:bin_width:3; %   [-1, 3] around cue onset, should work for both tasks
n_bins_cue = numel(histcounts([], bin_edges_cue));
%%  neuron class ANOVA
neuron_analysis_output_cue = struct([]);
% neuron_analysis_output_sac = struct([]);
for i = 1:numel(neuron_repo)
    if ~(neuron_tbl.task_id(i) == 1) %  Omit neuron if anti-saccade
        continue
    else
        tic
        n_class_ = numel(neuron_repo(i).class);
        if isempty(neuron_repo(i).class)
            continue
        end
        psth_cue_ = {neuron_repo(i).class.psth_cue};
%         psth_sac_ = {neuron_repo(i).class.psth_sac};
        for j = 1:n_bins_cue
            anova_y_cue_ = [];
            anova_label_cue_ = [];
%             anova_y_sac_ = [];
%             anova_label_sac_ = [];
            for m = 1:n_class_
                if m > 8
                    class_id_ = mod(m, 8);
                else
                    class_id_ = m;
                end
                anova_y_cue_ = [anova_y_cue_, psth_cue_{m}(:, j)'];
                anova_label_cue_ = [anova_label_cue_, class_id_ + zeros(size(psth_cue_{m}(:, j)'))];
%                 anova_y_sac_ = [anova_y_, psth_sac_{m}(:, j)'];
%                 anova_label_sac_ = [anova_label_, class_id_ + zeros(size(psth_sac_{m}(:, j)'))];
            end
            [~, out_tb_, out_stats_] = anova1(anova_y_cue_, anova_label_cue_, 'off');
            neuron_analysis_output_cue(i, j).tb = out_tb_;
            neuron_analysis_output_cue(i, j).stats = out_stats_;
%             [~, out_tb_, out_stats_] = anova1(anova_y_sac_, anova_label_sac_, 'off');
%             neuron_analysis_output_sac(i, j).tb = out_tb_;
%             neuron_analysis_output_sac(i, j).stats = out_stats_;
        end
        toc
    end
end
%%  Save neuron_analysis_output_cue
fname_ = 'neuron_analysis_output_cue_ps';
save(fullfile(project_dir, output_database, fname_), 'neuron_analysis_output_cue');
%%  Save neuron_analysis_output_sac
% fname_ = 'neuron_analysis_output_sac';
% save(fullfile(project_dir, output_database, fname_), 'neuron_analysis_output_sac');
%%  Compute PEV
neuron_pev_cue = zeros(size(neuron_analysis_output_cue));
% neuron_pev_sac = zeros(size(neuron_analysis_output_sac));
for i = 1:size(neuron_analysis_output_cue, 1)
%     i
    if neuron_tbl.task_id(i) == 1
        for j = 1:size(neuron_analysis_output_cue, 2)
            if isempty(neuron_analysis_output_cue(i, j).tb)
                        neuron_pev_cue(i, :) = NaN;        
        neuron_pev_cue(i, :) = NaN;
            else
                neuron_pev_cue(i, j) = zw_pev(neuron_analysis_output_cue(i, j).tb);
                if isnan(neuron_pev_cue(i, j))
                    i
                    j
                end
            end
        end
%         for j = 1:size(neuron_analysis_output_sac, 2)
%             neuron_pev_sac(i, j) = zw_pev(neuron_analysis_output_sac(i, j).tb);
%         end
    else
        neuron_pev_cue(i, :) = NaN;        
%         neuron_pev_sac(i, :) = NaN;
    end
end
%%  Save PEV
fname_ = 'neuron_pev_cue_ps';
save(fullfile(project_dir, output_database, fname_), 'neuron_pev_cue');
% fname_ = 'neuron_pev_sac';
% save(fullfile(project_dir, output_database, fname_), 'neuron_pev_sac');