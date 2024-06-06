clc; clear all; close all;
load('MSNG_rates.mat')
load('R1R2_rates.mat')
%%
[p, h, stats] = ranksum([MSNG_PFC_NS_cuerates; R1R2_PFC_NS_cuerates]', [MSNG_PFC_BS_cuerates; R1R2_PFC_BS_cuerates]')
[p, h, stats] = ranksum([MSNG_PPC_NS_cuerates; R1R2_PPC_NS_cuerates]', [MSNG_PPC_BS_cuerates; R1R2_PPC_BS_cuerates]')
[p, h, stats] = ranksum([MSNG_PFC_NS_fixrates; R1R2_PFC_NS_fixrates]', [MSNG_PFC_BS_fixrates; R1R2_PFC_BS_fixrates]')
[p, h, stats] = ranksum([MSNG_PPC_NS_fixrates; R1R2_PPC_NS_fixrates]', [MSNG_PPC_BS_fixrates; R1R2_PPC_BS_fixrates]')



%%
allNS_fix= [MSNG_PFC_NS_fixrates; R1R2_PFC_NS_fixrates; MSNG_PPC_NS_fixrates; R1R2_PPC_NS_fixrates];
allBS_fix= [MSNG_PFC_BS_fixrates; R1R2_PFC_BS_fixrates; MSNG_PPC_BS_fixrates; R1R2_PPC_BS_fixrates];
fprintf('Mean of NS population: %f\n', mean(allNS_fix))
fprintf('Mean of BS population: %f\n', mean(allBS_fix))

[p, h, stats] = ranksum(allNS_fix', allBS_fix')

[h, p, ci, stats] = ttest2(allNS_fix, allBS_fix)

%%
allNS_cue= [MSNG_PFC_NS_cuerates; R1R2_PFC_NS_cuerates; MSNG_PPC_NS_cuerates; R1R2_PPC_NS_cuerates];
allBS_cue= [MSNG_PFC_BS_cuerates; R1R2_PFC_BS_cuerates; MSNG_PPC_BS_cuerates; R1R2_PPC_BS_cuerates];
fprintf('Mean of NS population: %f\n', mean(allNS_cue))
fprintf('Mean of BS population: %f\n', mean(allBS_cue))

[p, h, stats] = ranksum(allNS_cue', allBS_cue')

[h, p, ci, stats] = ttest2(allNS_cue, allBS_cue)





