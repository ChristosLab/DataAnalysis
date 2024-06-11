%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 1}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 2}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 1}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 2}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 1}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 2}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 1}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 2}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 1}], 16:32, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 2}], 16:32, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 1}], 16:32, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 2}], 16:32, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 1}], 16:32, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 2}], 16:32, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 1}], 16:32, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 2}], 16:32, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 1}], 8:16, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 2}], 8:16, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 1}], 8:16, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 2}], 8:16, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 1}], 8:16, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 2}], 8:16, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 1}], 8:16, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 2}], 8:16, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 1}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 2}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 1}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 2}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 1}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 2}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 1}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 2}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 1}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 2}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 1}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 2}], 4:8, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
%%  Neuron PEV (nonsig)
pev_plot({neuron_pev_cue([neuron_tuning_cat_nonsig{1, [1,3,4], :, :}], :), ...
    neuron_pev_cue([neuron_tuning_cat_nonsig{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent nonsig', 'Adult nonsig'}, {[0.3, 0.85, 0.85], [0.85, 0.5, 0.3]})
title('Spiking PEV')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_spiking_tuning_draft'];
%%  Spiking firing rate raw (nonsig)
psth_plot({best_psth_raw([neuron_tuning_cat_nonsig{1, [1,3,4], :, :}], :), ...
    best_psth_raw([neuron_tuning_cat_nonsig{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent nonsig', 'Adult nonsig'}, {'b','r'})
title('Absolute Firing Rate')
xlabel('Time from cue onset (s)')
% ylim([0, 25])
fig_name = ['raw_spiking_firing_draft'];
%%  Neuron PEV by gamma power modulation (sig + nonsig)
pev_plot({neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 2}, neuron_mod_cat_nonsig{1, [1,3,4], 2}], :), ...
    neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 1}, neuron_mod_cat_nonsig{1, [1,3,4], 1}], :), ...
    neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 2}, neuron_mod_cat_nonsig{2, [1,3,4], 2}], :), ...
    neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 1}, neuron_mod_cat_nonsig{2, [1,3,4], 1}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Neuronal tuning')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_gamma_modulation'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Plot single trial bandpassed LFP
figure;
hold on
in = LFPData.class(4).ntr(1).LFP;
flp     = 32;
fhi     = 64;
fs      = 500;
wn      = [flp/(fs/2), fhi/(fs/2)];
[bb,aa] = butter(8,wn,'bandpass');
pl_freq = 60;
notch_w = 1/10;
stop_band = (pl_freq - 1 * (notch_w)):(notch_w):(pl_freq + 1 * (notch_w)); % Powerline noise removal.
in = in' - mean(in); %   Zero-mean
% plot(in)
filtered_ = ft_preproc_dftfilter(in,fs, stop_band);
filtered_2_ = filtfilt(bb,aa,filtered_);
% plot(filtered_);
plot(filtered_2_);
%%
%
%
%
%
%%  Spike field coherence
%%  Load spikes
fname_ = 'raster_repo_3_11_2021.mat';
load(fullfile(project_dir, output_database, fname_));
%%  Reduce repo size by emptying unwanted trials
for i = 1:numel(raster_repo)
    if ~ismember(i, [neuron_mod_cat{:}])
        raster_repo(i).class = [];
    end
end
%%  Load filed potentials
fname_ = 'lfp_repo_3_2_2021.mat';
load(fullfile(project_dir, output_database, fname_));
%%  Create segmented LFP repo
cue_dur = [-1, 4];
fs = 500;
lfp_segment_repo = struct('class', []);
tic
for i =1:numel(lfp_repo)
    lfp_segment_repo(i).class = struct();
    lfp_segment_repo(i).class.lfp_cue = zeros([0, 0, 0]);
end
for i = find(sum(mapping_mat_total, 2))'
    for j = 1:numel(lfp_repo(i).class)
        if isempty(lfp_repo(i).class(j).ntr)
            continue
        end
        %         if ~isfield(repo(i).class(j).ntr, 'Saccade_onT')
        %             continue
        %         end
        cue_sig = nan(numel(lfp_repo(i).class(j).ntr), round(diff(cue_dur)*fs));
        for k  =1:numel(lfp_repo(i).class(j).ntr)
            if isempty(lfp_repo(i).class(j).ntr(k).LFP)
                continue
            end
            %             if isempty(repo(i).class(j).ntr(k).Saccade_onT) || isnan(repo(i).class(j).ntr(k).Saccade_onT)
            %                 continue
            %             end
            cue_sample = floor(cue_dur*fs) + floor(lfp_repo(i).class(j).ntr(k).Cue_onT*fs);
            %             sac_sample = floor((sac_dur + repo(i).class(j).ntr(k).Saccade_onT)*fs);
            %
            %
            % Convert
            cue_sig(k, :) = lfp_repo(i).class(j).ntr(k).LFP((cue_sample(1) + 1):cue_sample(end))';
        end
        % Filtering
        cue_sig = zw_lfp_filt(cue_sig);
        lfp_segment_repo(i).class(j).lfp_cue = cue_sig;
    end
    toc
end
%%
fname_ = 'lfp_segment_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_segment_repo');
%%  Reduce repo size by emptying unwanted trials
for i = 1:numel(lfp_segment_repo)
    if ~ismember(i, [site_mod_cat{:}])
        lfp_segment_repo(i).class = [];
    end
end
%%  Spike-Field phase coherence gram to each neuron's preferred stimulus
target_frs = 2:2:128;
fs = 500;
down_sample = 20;
t_range_cue = [-1, 4];
cohgram_cue_repo = nan(numel(raster_repo), numel(target_frs), diff(t_range_cue)*fs/down_sample);
trial_counter = 0;
neuron_list_ = [neuron_tuning_cat{:, [1,3,4], :}];
mapping_mat_total = or(mapping_mat, mapping_mat_nonsig);
load(fullfile(project_dir, output_database, 'best_psth_upsampled.mat'), 'neuron_best_class');
%%
cue_dur = [-1, 4];
sf_cohegram = zeros(numel(raster_repo), numel(target_frs), diff(cue_dur)*1000/down_sample);
for i = neuron_list_
    tic
    lfp_id_ = find(mapping_mat_total(:, i));
    class_id_ = neuron_best_class(i);
    sig1 = raster_repo(i).class(class_id_).raster_cue_;
    sig2 = resample(lfp_segment_repo(lfp_id_).class(class_id_).lfp_cue', 1000, fs)';
    [coh_, ~] = zw_mtaper_cohegram(sig1, sig2, 1000, 1, 128);
    sf_cohegram(i, :, :) = downsample_row(squeeze(nanmean(coh_, 1)), down_sample);
    toc
end
%%
fname_ = 'sf_cohegram.mat';
save(fullfile(project_dir, output_database, fname_), 'sf_cohegram');
%%
figure
hold on
plot(target_frs, nanmean(nanmean(sf_cohegram([neuron_tuning_cat{1, :, 1}], :, 76:150), 3), 1))
plot(target_frs, nanmean(nanmean(sf_cohegram([neuron_tuning_cat{1, :, 2}], :, 76:150), 3), 1))
plot(target_frs, nanmean(nanmean(sf_cohegram([neuron_tuning_cat{2, :, 1}], :, 76:150), 3), 1))
plot(target_frs, nanmean(nanmean(sf_cohegram([neuron_tuning_cat{2, :, 2}], :, 76:150), 3), 1))
%%
figure
hold on
plot(target_frs, nanmean(nanmean(sf_cohegram([neuron_tuning_cat{1, :, 1}], :, 51:75), 3), 1))
plot(target_frs, nanmean(nanmean(sf_cohegram([neuron_tuning_cat{1, :, 2}], :, 51:75), 3), 1))
plot(target_frs, nanmean(nanmean(sf_cohegram([neuron_tuning_cat{2, :, 1}], :, 51:75), 3), 1))
plot(target_frs, nanmean(nanmean(sf_cohegram([neuron_tuning_cat{2, :, 2}], :, 51:75), 3), 1))
%%
figure
plot(squeeze(nanmean(nanmean(sf_cohegram([neuron_tuning_cat{1, :, :}], 8:16, :), 2), 1)))
hold on
plot(squeeze(nanmean(nanmean(sf_cohegram([neuron_tuning_cat{2, :, :}], 8:16, :), 2), 1)))
%%
fix_epoch   = 501:1000;
cue_epoch   = 1001:1500;
delay_epoch = 1501:3000;
sf_coherence_fix_epoch = zeros(numel(raster_repo), numel(target_frs));
sf_coherence_cue_epoch = zeros(numel(raster_repo), numel(target_frs));
sf_coherence_delay_epoch = zeros(numel(raster_repo), numel(target_frs));
for i = neuron_list_
    tic
    lfp_id_ = find(mapping_mat_total(:, i));
    class_id_ = neuron_best_class(i);
    sig1 = raster_repo(i).class(class_id_).raster_cue_;
    sig2 = resample(lfp_segment_repo(lfp_id_).class(class_id_).lfp_cue', 1000, fs)';
    coh_ = zw_mtaper_coherence(sig1(:, fix_epoch), sig2(:, fix_epoch), 1000, 2, 128);
    sf_coherence_fix_epoch(i, :) = nanmean(coh_, 1);
    coh_ = zw_mtaper_coherence(sig1(:, cue_epoch), sig2(:, cue_epoch), 1000, 2, 128);
    sf_coherence_cue_epoch(i, :) = nanmean(coh_, 1);
    coh_ = zw_mtaper_coherence(sig1(:, delay_epoch), sig2(:, delay_epoch), 1000, 2, 128);
    sf_coherence_delay_epoch(i, :) = nanmean(coh_, 1);
    toc
end
%%
figure
plot(target_frs, squeeze(nanmean(sf_coherence_fix_epoch([neuron_tuning_cat{1, :, 1}], :), 1)))
hold on
plot(target_frs, squeeze(nanmean(sf_coherence_fix_epoch([neuron_tuning_cat{1, :, 2}], :), 1)))
ylim([0, 0.45]);
%%
figure
plot(target_frs, squeeze(nanmean(sf_coherence_cue_epoch([neuron_tuning_cat{1, :, 1}], :), 1)))
hold on
plot(target_frs, squeeze(nanmean(sf_coherence_cue_epoch([neuron_tuning_cat{1, :, 2}], :), 1)))
ylim([0, 0.45]);
%%
figure
hold on
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch([neuron_tuning_cat{1, :, 1}], :), 1)))
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch([neuron_tuning_cat{1, :, 2}], :), 1)))
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch([neuron_tuning_cat{2, :, 1}], :), 1)))
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch([neuron_tuning_cat{2, :, 2}], :), 1)))
ylim([0, 0.45]);
%%
figure
hold on
for i = [neuron_tuning_cat{1, :, 2}]
    plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch(i, :), 1)))
    pause;
end
%%
test_lfp.label = {'lfp'};
test_lfp.trial = {1:1000};
test_lfp.time  = {0.001:0.001:1};
test_spike.timestamp = {random('Uniform', 0, 1, [1,10])};
test_spike.time = test_spike.timestamp;
test_spike.label = {'neuron1'};
test_spike.trial = {ones(1, 10)};
test_spike.trialtime = [0, 1];
test_out = ft_spiketriggeredspectrum(cfg, test_lfp,test_spike)
%%
fr_ = 20;
test_raster = zeros(size(sig2));
for i = 1:size(sig2, 1)
    test_raster(i, :) = zw_spike_time_to_raster(cumsum(exprnd(1/fr_, 1, size(sig2, 2)/1000*fr_), 2), [0, 5]);
end
test_out = zw_mtaper_coherence(test_raster(:, delay_epoch), sig2(:, delay_epoch), 1000, 2, 128);
figure; plot(squeeze(mean(test_out, 1)))
%%  Poisson process dummy
fix_epoch   = 501:1000;
cue_epoch   = 1001:1500;
delay_epoch = 1501:3000;
sf_coherence_fix_epoch_dummy = zeros(numel(raster_repo), numel(target_frs));
sf_coherence_cue_epoch_dummy = zeros(numel(raster_repo), numel(target_frs));
sf_coherence_delay_epoch_dummy = zeros(numel(raster_repo), numel(target_frs));
spike_cap = 200;
for i = neuron_list_
    tic
    lfp_id_ = find(mapping_mat_total(:, i));
    class_id_ = neuron_best_class(i);
    sig1 = raster_repo(i).class(class_id_).raster_cue_;
    sig2 = resample(lfp_segment_repo(lfp_id_).class(class_id_).lfp_cue', 1000, fs)';
    dum_n_spike_ = sum(sig1(:, fix_epoch), 2);
    dum_fr_      = dum_n_spike_/numel(fix_epoch)*1000;
    dum_raster_  = zw_spike_time_to_raster(cumsum(exprnd((1./dum_fr_).*ones(1, spike_cap), size(sig1, 1), spike_cap), 2), [0, 0.5]);
    coh_ = zw_mtaper_coherence(dum_raster_, sig2(:, fix_epoch), 1000, 2, 128);
    sf_coherence_fix_epoch_dummy(i, :) = nanmean(coh_, 1);
    dum_n_spike_ = sum(sig1(:, cue_epoch), 2);
    dum_fr_      = dum_n_spike_/numel(cue_epoch)*1000;
    dum_raster_  = zw_spike_time_to_raster(cumsum(exprnd((1./dum_fr_).*ones(1, spike_cap), size(sig1, 1), spike_cap), 2), [0, 0.5]);
    coh_ = zw_mtaper_coherence(dum_raster_, sig2(:, cue_epoch), 1000, 2, 128);
    sf_coherence_cue_epoch_dummy(i, :) = nanmean(coh_, 1);
    dum_n_spike_ = sum(sig1(:, delay_epoch), 2);
    dum_fr_      = dum_n_spike_/numel(delay_epoch)*1000;
    dum_raster_  = zw_spike_time_to_raster(cumsum(exprnd((1./dum_fr_).*ones(1, spike_cap), size(sig1, 1), spike_cap), 2), [0, 1.5]);
    coh_ = zw_mtaper_coherence(dum_raster_, sig2(:, delay_epoch), 1000, 2, 128);
    sf_coherence_delay_epoch_dummy(i, :) = nanmean(coh_, 1);
    toc
end
%%
figure
hold on
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch_dummy([neuron_tuning_cat{1, :, 1}], :), 1)))
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch_dummy([neuron_tuning_cat{1, :, 2}], :), 1)))
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch_dummy([neuron_tuning_cat{2, :, 1}], :), 1)))
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch_dummy([neuron_tuning_cat{2, :, 2}], :), 1)))
ylim([0, 0.45]);
%%  Periodic dummy
fix_epoch   = 501:1000;
cue_epoch   = 1001:1500;
delay_epoch = 1501:3000;
sf_coherence_fix_epoch_period   = zeros(numel(raster_repo), numel(target_frs));
sf_coherence_cue_epoch_period   = zeros(numel(raster_repo), numel(target_frs));
sf_coherence_delay_epoch_period = zeros(numel(raster_repo), numel(target_frs));
spike_cap = 200;
for i = neuron_list_
    tic
    lfp_id_ = find(mapping_mat_total(:, i));
    class_id_ = neuron_best_class(i);
    sig1 = raster_repo(i).class(class_id_).raster_cue_;
    sig2 = resample(lfp_segment_repo(lfp_id_).class(class_id_).lfp_cue', 1000, fs)';
    dum_n_spike_ = sum(sig1(:, fix_epoch), 2);
    dum_fr_      = dum_n_spike_/numel(fix_epoch)*1000;
    dum_raster_  = zw_spike_time_to_raster([1:spike_cap]./dum_fr_, [0, 0.5]);
    coh_ = zw_mtaper_coherence(dum_raster_, sig2(:, fix_epoch), 1000, 2, 128);
    sf_coherence_fix_epoch_period(i, :) = nanmean(coh_, 1);
    dum_n_spike_ = sum(sig1(:, cue_epoch), 2);
    dum_fr_      = dum_n_spike_/numel(cue_epoch)*1000;
    dum_raster_  = zw_spike_time_to_raster([1:spike_cap]./dum_fr_, [0, 0.5]);
    coh_ = zw_mtaper_coherence(dum_raster_, sig2(:, cue_epoch), 1000, 2, 128);
    sf_coherence_cue_epoch_period(i, :) = nanmean(coh_, 1);
    dum_n_spike_ = sum(sig1(:, delay_epoch), 2);
    dum_fr_      = dum_n_spike_/numel(delay_epoch)*1000;
    dum_raster_  = zw_spike_time_to_raster([1:spike_cap]./dum_fr_, [0, 1.5]);
    coh_ = zw_mtaper_coherence(dum_raster_, sig2(:, delay_epoch), 1000, 2, 128);
    sf_coherence_delay_epoch_period(i, :) = nanmean(coh_, 1);
    toc
end
%%
figure
hold on
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch_period([neuron_tuning_cat{1, :, 1}], :), 1)))
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch_period([neuron_tuning_cat{1, :, 2}], :), 1)))
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch_period([neuron_tuning_cat{2, :, 1}], :), 1)))
plot(target_frs, squeeze(nanmean(sf_coherence_delay_epoch_period([neuron_tuning_cat{2, :, 2}], :), 1)))
ylim([0, 0.45]);
%%
[c_, i_] = max(sf_coherence_delay_epoch_period([neuron_tuning_cat{1, :, 1}], :), [], 2);
figure
histogram(i_*2, 'BinEdges', [0:20:120], 'Normalization', 'probability')
%%
[c_, i_] = max(sf_coherence_delay_epoch_period([neuron_tuning_cat{1, :, 2}], :), [], 2);
figure
histogram(i_*2, 'BinEdges', [0:20:120], 'Normalization', 'probability')
%%
[c_, i_] = max(sf_coherence_delay_epoch([neuron_tuning_cat{1, :, :}], :), [], 2);
figure
histogram(i_*2, 'BinEdges', [0:10:120], 'Normalization', 'probability')
ylim([0, 0.25]);
%%
[c_, i_] = max(sf_coherence_delay_epoch([neuron_tuning_cat{1, :, 2}], :), [], 2);
figure
histogram(i_*2, 'BinEdges', [0:10:120], 'Normalization', 'probability')
ylim([0, 0.25]);
%%
[c_, i_] = max(sf_coherence_delay_epoch([neuron_tuning_cat{2, :, :}], :), [], 2);
figure
histogram(i_*2, 'BinEdges', [0:10:120], 'Normalization', 'probability')
ylim([0, 0.25]);
%%
[c_, i_] = max(sf_coherence_delay_epoch_dummy([neuron_tuning_cat{1, :, :}], :), [], 2);
figure
histogram(i_*2, 'BinEdges', [0:20:120], 'Normalization', 'probability')
%%
[c_, i_] = max(sf_coherence_delay_epoch_dummy([neuron_tuning_cat{2, :, :}], :), [], 2);
figure
histogram(i_*2, 'BinEdges', [0:20:120], 'Normalization', 'probability')
