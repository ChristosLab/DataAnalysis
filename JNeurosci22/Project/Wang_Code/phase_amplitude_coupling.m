zw_setpath;
addpath(fullfile(project_dir, code_lib, 'External\fieldtrip\'));
ft_defaults;
addpath(fullfile(project_dir, code_lib, 'External\PACmeg-master\functions\'));
%%
load(fullfile(project_dir, output_database, 'lfp_repo_complete.mat'));
load(fullfile(project_dir, output_database, 'lfp_tbl.mat'));
load(fullfile(project_dir, output_database, 'lfp_neuron_matching.mat'));
%%
t_1_t = intersect([t{1,[1,3,4]}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 2)));
t_1_nt = intersect([t{1,[1,3,4]}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 2)));
t_2_t = intersect([t{2,[1,3,4]}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 2)));
t_2_nt = intersect([t{2,[1,3,4]}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 2)));
%%
t_1_ = [t{1,[1,3,4]}];
t_2_ = [t{2,[1,3,4]}]';
%%
fs = 500;
% durs = [[1/fs, 0.5] + [-0.75:0.25:1.5]'; [0.5 + 1/fs, 2]];
% durs = [[1/fs, 0.5] + [-0.5:0.5:1.5]'; [0.5 + 1/fs, 2]];
durs = [0.5 + 1/fs, 2];
% durs = [0,5 + 1/fs, 1];
%
cfg_pac.Fs = fs;
cfg_pac.phase_freqs = 4:2:32;
cfg_pac.amp_freqs = 32:4:128;
% cfg_pac.amp_freqs = 60;
% cfg_pac.filt_order = 1;
cfg_pac.filt_order = 2;

% cfg_pac.amp_bandw_method = 'maxphase';
cfg_pac.amp_bandw_method = 'number';
cfg_pac.amp_bandw = 16;

cfg_pac.method = 'tort';
cfg_pac.durs = durs;

cfg_pac.n_surr = 500;
%%  Compute PAC including surrogate
target_lfp = [site_tuning_cat{:, [1,3,4], :, :}];
% zw_repo_pac(lfp_repo, cfg_pac, target_lfp);
pac_repo = zw_repo_pac(lfp_repo, cfg_pac, target_lfp);

%%
fname_ = 'pac_repo_5_6_21';
save(fullfile(project_dir, output_database, fname_), 'pac_repo', 'cfg_pac');
%%  Compute 1) surrogate p and 2) normalized PAC
for i = 1:numel(target_lfp)
    for j = 1:numel(pac_repo(target_lfp(i)).class)
        z_ = (pac_repo(target_lfp(i)).class(j).pac - pac_repo(target_lfp(i)).class(j).miu)./pac_repo(target_lfp(i)).class(j).delta;
        p_ = 1-normcdf(abs(z_));
        pac_repo(target_lfp(i)).class(j).normalized_pac = z_;
        pac_repo(target_lfp(i)).class(j).p = p_;
%         normalized_ = norminv(1 - pac_repo(target_lfp(i)).class(j).);
%         normalized_(isinf(normalized_)) = 8; % Z score cut-off
%         pac_repo(target_lfp(i)).class(j).normalized_pac = normalized_;
    end
end
%%  First-pass visual examination
for i = 1:numel(target_lfp)
    for j = 1:numel(pac_repo(target_lfp(i)).class)
        figure('Units', 'pixels', 'Position', [400, 400, 1000, 400]);
        subplot(1,2,1);
        imagesc(squeeze(pac_repo(target_lfp(i)).class(j).pac(1,:,:)), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs)
        colorbar
        subplot(1,2,2);
        imagesc(squeeze(pac_repo(target_lfp(i)).class(j).normalized_pac(1,:,:)), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [0, 5])

%         imagesc(squeeze(pac_repo(target_lfp(i)).class(j).p(1,:,:))<0.05, 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs)
        colorbar
        pause;
        close all
    end
end
%%
clm_n = 100;
clm_ = [[zeros(1, clm_n);zeros(1, clm_n); linspace(1, 0, clm_n)]'; [linspace(0, 1, clm_n); zeros(1, clm_n);zeros(1, clm_n)]'];
%%
z_pac_range =  [0.02, 0.62];
%%  Plot for each subgroup of sites: Mean normalized PAC
figure
imagesc(mean_pac(pac_repo(t_2_), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, z_pac_range)
comodulogram_setplot('Adult PAC', 'Surrogate-Normalized Modulation Index')
figure
imagesc(mean_pac(pac_repo(t_1_), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, z_pac_range)
comodulogram_setplot('Adolescent PAC', 'Surrogate-Normalized Modulation Index')
%%  DIFF
figure
imagesc(mean_pac(pac_repo(t_2_), 1, 1) - mean_pac(pac_repo(t_1_), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [-.35, .35])
colormap(clm_)
comodulogram_setplot('Adult PAC - Adolescent PAC', '\Delta Surrogate-Normalized Modulation Index')
%%  
max(mean_pac(pac_repo(t_2_), 1, 1) - mean_pac(pac_repo(t_1_), 1, 1), [], 'all')
%%  T TEST
figure
mat_ = ttest_pac(pac_repo(t_1_), pac_repo(t_2_), 1, 1);
imagesc(mat_, 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [0, .05])
colormap('gray')
comodulogram_setplot('Adult PAC ? Adolescent PAC', 'T-test uncorrected p value')
%%  Plot for each subgroup of sites: Mean raw PAC
figure
imagesc(mean_pac(pac_repo(t_1_), 1, 0), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs)
colorbar
%%  Normalized PAC young neuron no-tune/tune
figure
imagesc(mean_pac(pac_repo([site_tuning_cat{1, [1,3,4], :, 1}]), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, z_pac_range)
comodulogram_setplot('Adolescent PAC: Sites w/o tuned neurons', 'Surrogate-Normalized Modulation Index')
%
figure
imagesc(mean_pac(pac_repo([site_tuning_cat{1, [1,3,4], :, 2}]), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, z_pac_range)
comodulogram_setplot('Adolescent PAC: Sites w/ tuned neurons', 'Surrogate-Normalized Modulation Index')
%%
mat_ = mean_pac(pac_repo([site_tuning_cat{1, [1,3,4], :, 1}]), 1, 1);
m_ = max(mat_, [], 'All')
[rowsOfMaxes colsOfMaxes] = find(mat_ == m_)
mat_ = mean_pac(pac_repo([site_tuning_cat{1, [1,3,4], :, 2}]), 1, 1);
m_ = max(mat_, [], 'All')
[rowsOfMaxes colsOfMaxes] = find(mat_ == m_)
%%  Normalized PAC adult neuron no-tune/tune
figure
imagesc(mean_pac(pac_repo([site_tuning_cat{2, [1,3,4], :, 1}]), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, z_pac_range)
comodulogram_setplot('Adult PAC: Sites w/o tuned neurons', 'Surrogate-Normalized Modulation Index')
%
figure
imagesc(mean_pac(pac_repo([site_tuning_cat{2, [1,3,4], :, 2}]), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs,  z_pac_range)
comodulogram_setplot('Adult PAC: Sites w/ tuned neurons', 'Surrogate-Normalized Modulation Index')
%%
mat_ = mean_pac(pac_repo([site_tuning_cat{2, [1,3,4], :, 1}]), 1, 1);
m_ = max(mat_, [], 'All')
[rowsOfMaxes colsOfMaxes] = find(mat_ == m_)
mat_ = mean_pac(pac_repo([site_tuning_cat{2, [1,3,4], :, 2}]), 1, 1);
m_ = max(mat_, [], 'All')
[rowsOfMaxes colsOfMaxes] = find(mat_ == m_)
%%
mat_ = sd_pac(pac_repo([site_tuning_cat{2, [1,3,4], :, 2}]), 1, 1);
%%  t-test adult neuron tune vs no-tune
mat_ = ttest_pac(pac_repo([site_tuning_cat{2, [1,3,4], :, 1}]), pac_repo([site_tuning_cat{2, [1,3,4], :, 2}]), 1, 1);
figure
imagesc(mat_, 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs)
colorbar
%%  t-test adult neuron tune vs adolescent neuron tune
mat_ = ttest_pac(pac_repo([site_tuning_cat{2, [1,3,4], :, 2}]), pac_repo([site_tuning_cat{1, [1,3,4], :, 2}]), 1, 1);
figure
imagesc(mat_, 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs)
colorbar
%%  t-test adult neuron all vs adolescent neuron all
mat_ = ttest_pac(pac_repo(t_1_), pac_repo(t_2_), 1, 1);
figure
imagesc(mat_ < 0.05, 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs)
colorbar
%%  Normalized PAC young gamma no-tune/tune
figure
imagesc(mean_pac(pac_repo([site_tuning_cat{1, [1,3,4], 1, :}]), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [0.7, 1.15])
colorbar
%
figure
imagesc(mean_pac(pac_repo([site_tuning_cat{1, [1,3,4], 2, :}]), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [0.7, 1.15])
colorbar
%%  Normalized PAC adult gamma no-tune/tune
figure
imagesc(mean_pac(pac_repo([site_tuning_cat{2, [1,3,4], 1, :}]), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [0.7, 1.15])
colorbar
%
figure
imagesc(mean_pac(pac_repo([site_tuning_cat{2, [1,3,4], 2, :}]), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [0.7, 1.15])
colorbar
%%  Normalized pac adult_tune - young_tune
figure
imagesc(mean_pac(pac_repo(t_2_t), 1, 1) - mean_pac(pac_repo(t_1_t), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [-.3, .3])
colormap(clm_)
colorbar
%%  Normalized pac adult_all - young_all
figure
imagesc(mean_pac(pac_repo(t_2_), 1, 1) - mean_pac(pac_repo(t_1_), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [-.3, .3])
colormap(clm_)
colorbar

%%  Raw pac percentage adult_all - young_all
figure
imagesc((mean_pac(pac_repo(t_2_), 1, 0) - mean_pac(pac_repo(t_1_), 1, 0))./mean_pac(pac_repo(t_1_), 1, 0), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs)
colorbar
%%  Raw pac percentage adult_tune - young_tune
figure
imagesc((mean_pac(pac_repo(t_2_nt), 1, 0) - mean_pac(pac_repo(t_1_nt), 1, 0))./mean_pac(pac_repo(t_1_), 1, 0), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs)
colorbar
%%  Raw pac percentage adult_del_mod - young_del_mod
figure
imagesc((mean_pac(pac_repo(intersect(t_2_, find(tfr_delay_modulation_flag(:, 3, 1)))), 1, 0) - mean_pac(pac_repo(intersect(t_1_, find(tfr_delay_modulation_flag(:, 3, 1)))), 1, 0))./mean_pac(pac_repo(intersect(t_1_, find(tfr_delay_modulation_flag(:, 3, 1)))), 1, 0), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [-0.8, 0.8])
colormap(clm_)
colorbar
%%  Raw pac percentage adult_del_no_mod - young_del_no_mod
figure
imagesc((mean_pac(pac_repo(intersect(t_2_, find(~tfr_delay_modulation_flag(:, 3, 1)))), 1, 0) - mean_pac(pac_repo(intersect(t_1_, find(~tfr_delay_modulation_flag(:, 3, 1)))), 1, 0))./mean_pac(pac_repo(intersect(t_1_, find(~tfr_delay_modulation_flag(:, 3, 1)))), 1, 0), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [-0.8, 0.8])
colormap(clm_)
colorbar
%%  Normalized pac adult_del_mod - young_del_mod
figure
imagesc(mean_pac(pac_repo(intersect(t_2_, find(tfr_delay_modulation_flag(:, 3, 1)))), 1, 1) - mean_pac(pac_repo(intersect(t_1_, find(tfr_delay_modulation_flag(:, 3, 1)))), 1, 1), 'XData', cfg_pac.phase_freqs, 'YData', cfg_pac.amp_freqs, [-0.4, 0.4])
colormap(clm_)
colorbar
%%
%%
zw_anova1_from_cell()
%%
pac_for_plot = zeros(0,size(cfg_pac.durs, 1),numel(cfg_pac.amp_freqs),numel(cfg_pac.phase_freqs));
for i = 1:numel(pac_repo)
    if ~any([t_1_, t_2_] == i)
        continue
    end
    pac_ = zeros(0, 10, 48, 13);
    for k = 2:numel(pac_repo(i).class)
        pac_ = [pac_; pac_repo(i).class(k).pac];
    end
        pac_for_plot(i, :, :, :) = nanmean(pac_, 1);
end
%%
function out_mean = mean_pac(in_repo, window, normalized)
mean_ = pac_session_mean(in_repo, window, normalized);
out_mean = squeeze(mean(mean_, 1));
end
%%
function out_sd = sd_pac(in_repo, window, normalized)
mean_ = pac_session_mean(in_repo, window, normalized);
out_sd = squeeze(std(mean_, 1));
end
%%
function out_p = ttest_pac(in_repo_1, in_repo_2, window, normalized)
mean_1 = pac_session_mean(in_repo_1, window, normalized);
mean_2 = pac_session_mean(in_repo_2, window, normalized);
out_p = zeros(size(mean_1, 2), size(mean_1, 3));
for i = 1:size(mean_1, 2)
    for j =1:size(mean_1, 3)
        [~, p_] = ttest2(mean_1(:, i, j), mean_2(:, i, j));
        out_p(i, j) = p_;
    end
end
end 
%%
function out = pac_session_mean(in_repo, window, normalized)
out = zeros([numel(in_repo), size(squeeze(in_repo(1).class(1).normalized_pac(window, :, :)))]);
for i = 1:numel(in_repo)
    for_plot = zeros([0, size(squeeze(in_repo(1).class(1).normalized_pac(window, :, :)))]);
    for j = 1:numel(in_repo(i).class)
        if normalized
            for_plot = [for_plot; in_repo(i).class(j).normalized_pac(window, :, :)];
        else
            for_plot = [for_plot; in_repo(i).class(j).pac(window, :, :)];
        end
        out(i, :, :) = mean(for_plot, 1);
    end
end
end
%%
function comodulogram_setplot(title_str, colorbar_str)
    xlabel('Phase Frequency (Hz)')
    ylabel('Amplitude Frequency (Hz)')
    set(gca, 'YDir', 'normal');
    h = colorbar;
    ylabel(h, colorbar_str);
    title(title_str);
    set(gca, 'FontSize', 14);
    set(gca, 'FontWeight', 'bold');
    fig_name = matlab.lang.makeValidName(strcat(title_str, string(datetime(now,'ConvertFrom','datenum'))));
    print(gcf,  fullfile('D:\Database\Zhengyang_Wang\Projects\LFP_project_v0.2\Wang_Figure', fig_name), '-dpng', '-r400');
end