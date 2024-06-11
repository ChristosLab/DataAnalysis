figure(); plot(mean(squeeze(nanmean(temp_cwt(t11, 12:17, :)./temp_baseline(t11, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t21, 12:17, :)./temp_baseline(t21, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt(t11, 32:64, :)./temp_baseline(t11, 32:64), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t21, 32:64, :)./temp_baseline(t21, 32:64), 1)), 1));%ylim([6, 8]*10^-3)
figure(); plot(mean(squeeze(nanmean(temp_cwt(t11, 8:16, :)./temp_baseline(t11, 8:16), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt(t21, 8:16, :)./temp_baseline(t21, 8:16), 1)), 1));%ylim([6, 8]*10^-3)
%%
figure()
histogram(nanmean(nanmean(temp_cwt(t11, 4:8, 51:75)./temp_baseline(t11, 4:8), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.1:2)
hold on

histogram(nanmean(nanmean(temp_cwt(t21, 4:8, 51:75)./temp_baseline(t21, 4:8), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.1:2)
%%
figure()
histogram(nanmean(nanmean(temp_cwt(t11, 8:16, 51:75)./temp_baseline(t11, 8:16), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
hold on

histogram(nanmean(nanmean(temp_cwt(t21, 8:16, 51:75)./temp_baseline(t21, 8:16), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
legend({'Adolescent', 'Adult'})
%%
figure()
histogram(nanmean(nanmean(temp_cwt(t11, 32:64, 51:75)./temp_baseline(t11, 32:64), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
hold on

histogram(nanmean(nanmean(temp_cwt(t21, 32:64, 51:75)./temp_baseline(t21, 32:64), 2), 3), 'Normalization', 'probability', 'BinEdges', 0.5:0.05:2)
legend({'Adolescent', 'Adult'})