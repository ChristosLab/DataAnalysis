figure(); plot(mean(squeeze(nanmean(temp_cwt(t13, 12:17, :)./temp_baseline(t13, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t23, 12:17, :)./temp_baseline(t23, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt(t13, 32:64, :)./temp_baseline(t13, 32:64), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t23, 32:64, :)./temp_baseline(t23, 32:64), 1)), 1));%ylim([6, 8]*10^-3)
figure(); plot(mean(squeeze(nanmean(temp_cwt(t13, 8:16, :)./temp_baseline(t13, 8:16), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t23, 8:16, :)./temp_baseline(t23, 8:16), 1)), 1));%ylim([6, 8]*10^-3)
%%
figure();hist(nanmean(nanmean(temp_cwt(t23, 8:16, 51:75)./temp_baseline(t23, 8:16), 2), 3))
%%
figure();hist(nanmean(nanmean(temp_cwt(t23,32:64, 51:75)./temp_baseline(t23, 32:64), 2), 3), 0.5:0.1:2)
%%
figure();hist(nanmean(nanmean(temp_cwt(t23,16:32, 51:75)./temp_baseline(t23, 16:32), 2), 3), 0.5:0.1:2)
%%
figure();hist(nanmean(nanmean(temp_cwt(t13, 8:16, 51:75)./temp_baseline(t13, 8:16), 2), 3))
%%
figure();hist(nanmean(nanmean(temp_cwt(t13,32:64, 51:75)./temp_baseline(t13, 32:64), 2), 3), 0.5:0.1:2)
%%
figure();hist(nanmean(nanmean(temp_cwt(t13,16:32, 51:75)./temp_baseline(t13, 16:32), 2), 3), 0.5:0.1:2)
%%
figure();hist(nanmean(nanmean(temp_cwt(t13,4:8, 51:75)./temp_baseline(t13, 4:8), 2), 3), 0.5:0.1:2)
%%
figure()
histogram(nanmean(nanmean(temp_cwt(t13, 4:8, 51:75)./temp_baseline(t13, 4:8), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.1:2)
hold on

histogram(nanmean(nanmean(temp_cwt(t23, 4:8, 51:75)./temp_baseline(t23, 4:8), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.1:2)
%%
figure()
histogram(nanmean(nanmean(temp_cwt(t13, 8:16, 51:75)./temp_baseline(t13, 8:16), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
hold on

histogram(nanmean(nanmean(temp_cwt(t23, 8:16, 51:75)./temp_baseline(t23, 8:16), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
legend({'Adolescent', 'Adult'})
%%
figure()
histogram(nanmean(nanmean(temp_cwt(t13, 32:64, 51:75)./temp_baseline(t13, 32:64), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
hold on

histogram(nanmean(nanmean(temp_cwt(t23, 32:64, 51:75)./temp_baseline(t23, 32:64), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
legend({'Adolescent', 'Adult'})