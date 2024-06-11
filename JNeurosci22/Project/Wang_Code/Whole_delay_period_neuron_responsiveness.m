%%  Whole delay period neuron responsiveness
fname_ = 'complete_neuron_repo_w_psth.mat';
load(fullfile(project_dir, output_database, fname_));
%%
% n_bins_lfp = (diff(sac_dur)*fs + 1)/down_sample;
neuron_epoch_fr_cue = struct([]);
% lfp_analysis_output_sac = struct([]);
tic
for i = 1:numel(neuron_repo)
    if ~(neuron_tbl.task_id(i) == 1) %  Omit neuron if anti-saccade
        neuron_epoch_fr_cue(i).mean = [];
        neuron_epoch_fr_cue(i).p = [];
        continue
    end
    if isempty(neuron_repo(i).class)
        neuron_epoch_fr_cue(i).mean = [];
        neuron_epoch_fr_cue(i).p = [];
        continue
    end
    n_class_ = numel(neuron_repo(i).class);
    cue_neuron_epoch_ = {neuron_repo(i).class.psth_cue};
    fr_ = [];
    for m = 1:n_class_
        if m > 8
            class_id_ = mod(m, 8);
        else
            class_id_ = m;
        end
        fr_ = [fr_, nanmean(cue_neuron_epoch_{m}(:, 31:60), 2)'];
    end
    [h_, p_] = ttest(fr_, 0);
            neuron_epoch_fr_cue(i).mean = nanmean(fr_);
            neuron_epoch_fr_cue(i).p = p_;
    toc
end
%%
fname_ = 'neuron_epoch_delay_responsiveness.mat';
save(fullfile(project_dir, output_database, fname_), 'neuron_epoch_fr_cue');
