t1 = [t13;t14]; %   Adolescent LFP sites
t2 = [t23;t24]; %   Adult LFP sites
%%
figure()
plot(nanmean(average_neuron_pev_by_site([t13;t14], 40:50), 2), nanmean(lfp_pev_cue([t13;t14], 100:125, 3), 2), '.')
hold on
plot(nanmean(average_neuron_pev_by_site([t23;t24], 40:50), 2), nanmean(lfp_pev_cue([t23;t24], 100:125, 3), 2), '.')
%%
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
%% Cue: LFP tuning vs neuron tuning 
xl = [-.4, 1.2; -.3, .5; -.4, .4; -.4, .8];
for lfp_band = 1:4
neuron_time = 20:30;
lfp_time = (neuron_time(1)*2.5):(neuron_time(2)*2.5);
LFP_PEV = nanmean(lfp_pev_cue([t1;t2], lfp_time, lfp_band), 2);
stage = [1*ones(size(t1)); 2*ones(size(t2))];
neuron_PEV = nanmean(average_neuron_pev_by_site([t1;t2], neuron_time), 2);
glm_tb = table(LFP_PEV, neuron_PEV, stage);
glm_tb.stage = categorical(glm_tb.stage);
mspec = 'neuron_PEV~LFP_PEV*stage';
mdl = fitlm(glm_tb, mspec);
figure()
gscatter(LFP_PEV, neuron_PEV,stage, 'br','xo')
w = linspace(min(LFP_PEV),max(LFP_PEV));
line(w,feval(mdl,w,'1'),'Color','b','LineWidth',1.5)
line(w,feval(mdl,w,'2'),'Color','r','LineWidth',1.5)
title([title_st{lfp_band}, '  ANCOVA Adjusted {\it R^2} = ', num2str(mdl.Rsquared.Adjusted, '%.3f')])
legend({'Adolescent', 'Adult'})
xlabel('LFP PEV');
ylabel('Neuron PEV');
% xlim(xl(lfp_band, :));
set_plot_12_4();
fig_name = ['lfp_tune_neuron_tune_cue_band_', num2str(lfp_band)];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%% Delay: LFP tuning vs neuron tuning 
xl = [-.7, .7; -.65, .8; -.6, 1.2; -.8, 1.2];
for lfp_band = 1:4
neuron_time = 30:60;
lfp_time = (neuron_time(1)*2.5):(neuron_time(2)*2.5);
LFP_PEV = nanmean(lfp_pev_cue([t1;t2], lfp_time, lfp_band), 2);
stage = [1*ones(size(t1)); 2*ones(size(t2))];
neuron_PEV = nanmean(average_neuron_pev_by_site([t1;t2], neuron_time), 2);
glm_tb = table(LFP_PEV, neuron_PEV, stage);
glm_tb.stage = categorical(glm_tb.stage);
mspec = 'neuron_PEV~LFP_PEV*stage';
mdl = fitlm(glm_tb, mspec);
figure()
gscatter(LFP_PEV, neuron_PEV,stage, 'br','xo')
w = linspace(min(LFP_PEV),max(LFP_PEV));
line(w,feval(mdl,w,'1'),'Color','b','LineWidth',1.5)
line(w,feval(mdl,w,'2'),'Color','r','LineWidth',1.5)
title([title_st{lfp_band}, '  ANCOVA Adjusted {\it R^2} = ', num2str(mdl.Rsquared.Adjusted, '%.3f')])
legend({'Adolescent', 'Adult'})
xlabel('LFP PEV');
ylabel('Neuron PEV');
% xlim(xl(lfp_band, :));
set_plot_12_4();
fig_name = ['lfp_tune_neuron_tune_delay_band_', num2str(lfp_band)];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%% Cue: LFP modulation vs neuron tuning 
xl = [-.4, 1.2; -.3, .5; -.4, .4; -.4, .8];
for lfp_band = 1:4
neuron_time = 20:30;
lfp_time = (neuron_time(1)*2.5):(neuron_time(2)*2.5);
LFP_pow = squeeze(nanmean(nanmean(temp_cwt([t1;t2], target_frs(lfp_band, 1):target_frs(lfp_band, 2), lfp_time), 3), 2)./nanmean(temp_baseline([t1;t2], target_frs(lfp_band, 1):target_frs(lfp_band, 2)), 2));
LFP_pow = LFP_pow - 1;
stage = [1*ones(size(t1)); 2*ones(size(t2))];
neuron_PEV = nanmean(average_neuron_pev_by_site([t1;t2], neuron_time), 2);
glm_tb = table(LFP_pow, neuron_PEV, stage);
glm_tb.stage = categorical(glm_tb.stage);
mspec = 'neuron_PEV~LFP_pow*stage';
mdl = fitlm(glm_tb, mspec);
figure()
gscatter(LFP_pow, neuron_PEV,stage, 'br','xo')
w = linspace(min(LFP_pow),max(LFP_pow));
line(w,feval(mdl,w,'1'),'Color','b','LineWidth',1.5)
line(w,feval(mdl,w,'2'),'Color','r','LineWidth',1.5)
title([title_st{lfp_band}, '  ANCOVA Adjusted {\it R^2} = ', num2str(mdl.Rsquared.Adjusted, '%.3f')])
legend({'Adolescent', 'Adult'})
xlabel('LFP power modulation');
ylabel('Neuron PEV');
xlim(xl(lfp_band, :));
set_plot_12_4();
fig_name = ['lfp_mod_neuron_tune_cue_band_', num2str(lfp_band)];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%% Delay: LFP modulation vs neuron tuning 
xl = [-.7, .7; -.65, .8; -.6, 1.2; -.8, 1.2];
for lfp_band = 1:4
neuron_time = 30:60;
lfp_time = (neuron_time(1)*2.5):(neuron_time(2)*2.5);
LFP_pow = squeeze(nanmean(nanmean(temp_cwt([t1;t2], target_frs(lfp_band, 1):target_frs(lfp_band, 2), lfp_time), 3), 2)./nanmean(temp_baseline([t1;t2], target_frs(lfp_band, 1):target_frs(lfp_band, 2)), 2));
LFP_pow = LFP_pow - 1;
stage = [1*ones(size(t1)); 2*ones(size(t2))];
neuron_PEV = nanmean(average_neuron_pev_by_site([t1;t2], neuron_time), 2);
glm_tb = table(LFP_pow, neuron_PEV, stage);
glm_tb.stage = categorical(glm_tb.stage);
mspec = 'neuron_PEV~LFP_pow*stage';
mdl = fitlm(glm_tb, mspec);
figure()
gscatter(LFP_pow, neuron_PEV,stage, 'br','xo')
w = linspace(min(LFP_pow),max(LFP_pow));
line(w,feval(mdl,w,'1'),'Color','b','LineWidth',1.5)
line(w,feval(mdl,w,'2'),'Color','r','LineWidth',1.5)
title([title_st{lfp_band}, '  ANCOVA Adjusted {\it R^2} = ', num2str(mdl.Rsquared.Adjusted, '%.3f')])
legend({'Adolescent', 'Adult'})
xlabel('LFP power modulation');
ylabel('Neuron PEV');
xlim(xl(lfp_band, :));
set_plot_12_4();
fig_name = ['lfp_mod_neuron_tune_delay_band_', num2str(lfp_band)];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%%
neuron_time = 30:60;
lfp_time = (neuron_time(1)*2.5):(neuron_time(2)*2.5);
lfp_band = 4;
LFP_pow = squeeze(nanmean(nanmean(temp_cwt([t1;t2], target_frs(lfp_band, 1):target_frs(lfp_band, 2), lfp_time), 3), 2)./nanmean(temp_baseline([t1;t2], target_frs(lfp_band, 1):target_frs(lfp_band, 2)), 2));
LFP_pow = LFP_pow - 1;
stage = [1*ones(size(t1)); 2*ones(size(t2))];
LFP_PEV = nanmean(lfp_pev_cue([t1;t2], lfp_time, lfp_band), 2);
glm_tb = table(LFP_pow, LFP_PEV, stage);
glm_tb.stage = categorical(glm_tb.stage);
mspec = 'LFP_PEV~LFP_pow*stage';
mdl = fitlm(glm_tb, mspec)
figure()
gscatter(LFP_pow, LFP_PEV,stage, 'br','xo')
legend({'Adolescent', 'Adult'})
% line(w,feval(fit,w,'70'),'Color','b','LineWidth',2)
% line(w,feval(fit,w,'76'),'Color','g','LineWidth',2)
% line(w,feval(fit,w,'82'),'Color','r','LineWidth',2)