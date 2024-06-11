[temp_wt, temp_f] = inspect_plot(LFP_repo, 10);
%%
fun2 = @(x) max(LFP_repo(x).LFP);
temp_range = arrayfun(fun2, 1:numel(LFP_repo), 'UniformOutput', false);
%%
function [wt, f] = inspect_plot(LFP_repo, ind)
prep_dur = 1;
n_sample = prep_dur*LFP_repo(ind).Sampling_Rate;
figure();
[wt, f] = cwt(LFP_repo(ind).LFP - mean(LFP_repo(ind).LFP), LFP_repo(ind).Sampling_Rate);
figure();
plot((0:length(LFP_repo(ind).LFP)-1)./LFP_repo(ind).Sampling_Rate, LFP_repo(ind).LFP);
end