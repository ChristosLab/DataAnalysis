%%  Create a neuron repo with all available neurons
%   'neuron_tbl.mat' was created with script 'create_neuron_tbl.m', which
%   contains identity information of all 5538 neurons available.
%   'neuron_repo_w_psth.mat' and 'needed_neuron.mat' was created in script 
%   'run_repo_cwt_9_14.m', which contains only
%   neurons assoicated with the first batch of LFP files
%%
clear
zw_setpath;
%%  Load neuron table
fname_ = 'neuron_tbl.mat';
load(fullfile(project_dir, output_database, fname_));
fname_ = 'needed_neuron.mat';
load(fullfile(project_dir, output_database, fname_));
fname_ = 'neuron_repo_w_psth.mat';
load(fullfile(project_dir, output_database, fname_));
%%  Save neuron data in one repo
% file_names_ = strcat(neuron_tbl.Filename(find(needed_neuron)), '_', num2str(neuron_tbl.Neuron(find(needed_neuron))));
file_names_ = strcat(neuron_tbl.Filename, '_', num2str(neuron_tbl.Neuron));
for i = 1:numel(file_names_)
    if ~needed_neuron(i)
        try
            load(fullfile(neuron_database, file_names_{i}))
            neuron_repo(i) = MatData;
            clear MatData
        catch
            disp(['error loading neuron file: ', file_names_{i}, newline])
            needed_neuron(i) = 0;
        end
        i
    end
end
%%
fname_ = 'compelte_neuron_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'neuron_repo');
%%  Check for non-exising neuron files
m = 0;
for i =1:numel(neuron_repo)
    if numel(neuron_repo(i).class) < 1
        m = m + 1
    end
end
%%
fname_ = 'compelte_neuron_repo.mat';
load(fullfile(project_dir, output_database, fname_));
%%  compute psth either cue or saccade aligned
no_sac = [];
bin_width = 0.05;  % 50 milliseconds bin
bin_edges_cue = -1:bin_width:3; %   [-1, 3] around cue onset, should work for both tasks
bin_edges_sac = -3:bin_width:1; %   [-3, 1] around sccade onset
t_bins_cue = bin_edges_cue(1:(end-1)) + bin_width/2;
t_bins_sac = bin_edges_sac(1:(end-1)) + bin_width/2;
n_bins_cue = numel(histcounts([], bin_edges_cue));
n_bins_sac = numel(histcounts([], bin_edges_sac));
for i = find(~needed_neuron)'
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
            psth_cue_(m, :) = histcounts(TS_cue_, bin_edges_cue)/bin_width;
            if isfield(neuron_repo(i).class(j).ntr(m), 'Saccade_onT')
                if ~isempty(neuron_repo(i).class(j).ntr(m).Saccade_onT)
                    TS_sac_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Saccade_onT;
                    psth_sac_(m, :) = histcounts(TS_sac_, bin_edges_sac)/bin_width;
                else
                    psth_sac_(m, :) = NaN;
                    no_sac = [no_sac, i];
                end
            else
                psth_sac_(m, :) = NaN;
                no_sac = [no_sac, i];
            end
        end
        neuron_repo(i).class(j).psth_cue = psth_cue_;
        neuron_repo(i).class(j).psth_sac = psth_sac_;
    end
end
%%
no_sac_u = size(unique(no_sac));
no_sac_u = unique(no_sac);
no_sac_u_1 = zeros(size(needed_neuron));
for i = no_sac_u
no_sac_u_1(i) = 1;
end
%%
fname_ = 'complete_neuron_repo_w_psth.mat';
save(fullfile(project_dir, output_database, fname_), 'neuron_repo');