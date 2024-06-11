fname_ = 'upsample_neuron_repo.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_repo');
%%
raster_repo = struct;
t_range_cue = [-1, 4];
t_range_sac = [-2, 3];
step = 0.001;
bin_edges_cue = t_range_cue(1):step:t_range_cue(2);
bin_edges_sac = t_range_sac(1):step:t_range_sac(2);
n_bins_cue = numel(histcounts([], bin_edges_cue));
n_bins_sac = numel(histcounts([], bin_edges_sac));
neuron_list_ = [neuron_tuning_cat_nonsig{:, [1,3,4], :}, neuron_tuning_cat{:, [1,3,4], :}];
% for i = neuron_list_ % All LFP matched sig and non-sig neurons
for i = 1:numel(neuron_repo)
    if numel(neuron_repo(i).class) < 1
        continue
    end
    i
    for j = 1:numel(neuron_repo(i).class)
        if isempty(neuron_repo(i).class(j).ntr)
            continue
        end
        raster_cue_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_cue, 'logical');
        raster_sac_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_sac, 'logical');
        for m = 1:numel(neuron_repo(i).class(j).ntr)
            TS_cue_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Cue_onT;
            raster_cue_(m, :) = zw_spike_time_to_raster(TS_cue_, t_range_cue);
            if isfield(neuron_repo(i).class(j).ntr(m), 'Saccade_onT')
                if ~isempty(neuron_repo(i).class(j).ntr(m).Saccade_onT)
                    TS_sac_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Saccade_onT;
                    raster_sac_(m, :) = zw_spike_time_to_raster(TS_sac_, t_range_sac);
                end
            end
        end
        raster_repo(i).class(j).raster_cue_ = raster_cue_;
        raster_repo(i).class(j).raster_sac_ = raster_sac_;
    end
end
%%
fname_ = 'raster_repo_3_11_2021.mat';
save(fullfile(project_dir, output_database, fname_), 'raster_repo');