bin_width = 0.05;  % 50 milliseconds bin
step = 0.02; %  20 millisecond step matching that of SPD
t_range_cue = [-1, 3]; %   [-1, 3] around cue onset, should work for both tasks
t_range_sac = [-3, 1];
bin_edges_cue = -1:step:3;
bin_edges_sac = -3:step:1;
n_bins_cue = numel(histcounts([], bin_edges_cue));
n_bins_sac = numel(histcounts([], bin_edges_sac));
for i = [neuron_tuning_cat{:, [1,3,4], :}]
    if numel(neuron_repo(i).class) < 1
        continue
    end
    i
    for j = 1:numel(neuron_repo(i).class)
        if isempty(neuron_repo(i).class(j).ntr)
            continue
        end
        psth_cue_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_cue);
        psth_sac_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_sac);
        for m = 1:numel(neuron_repo(i).class(j).ntr)
            TS_cue_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Cue_onT;
            psth_cue_(m, :) = zw_spike_time_to_psth(TS_cue_, bin_width, step, t_range_cue);
            if isfield(neuron_repo(i).class(j).ntr(m), 'Saccade_onT')
                if ~isempty(neuron_repo(i).class(j).ntr(m).Saccade_onT)
                    TS_sac_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Saccade_onT;
                    psth_sac_(m, :) = zw_spike_time_to_psth(TS_sac_, bin_width, step, t_range_sac);
                else
                    psth_sac_(m, :) = NaN;
                end
            else
                psth_sac_(m, :) = NaN;
            end
        end
        neuron_repo(i).class(j).psth_cue_upsampled = psth_cue_;
        neuron_repo(i).class(j).psth_sac_upsampled = psth_sac_;
    end
end
%%
fname_ = 'upsample_neuron_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'neuron_repo');