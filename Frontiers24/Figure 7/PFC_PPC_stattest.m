clear all; close all; clc;

load("PFC_roc.mat")

load("PPC_roc.mat")

%%
disp('PFC NS delay1 vs PPC NS delay1')
[h,p,ci,stats] = ttest2(PFC_delay_NS, PPC_delay_NS)
meanDifference = mean(PFC_delay_NS) - mean(PPC_delay_NS);
pooledStdDev = sqrt(((std(PFC_delay_NS)^2 + std(PPC_delay_NS)^2) / 2));
cohenD = meanDifference / pooledStdDev

disp('PFC BS delay1 vs PPC BS delay1')
[h,p,ci,stats] = ttest2(PFC_delay_BS, PPC_delay_BS)
meanDifference = mean(PFC_delay_BS) - mean(PPC_delay_BS);
pooledStdDev = sqrt(((std(PFC_delay_BS)^2 + std(PPC_delay_BS)^2) / 2));
cohenD = meanDifference / pooledStdDev

















