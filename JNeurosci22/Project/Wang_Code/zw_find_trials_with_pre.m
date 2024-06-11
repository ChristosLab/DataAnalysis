function [pre, post] = zw_find_trials_with_pre(in_trial_num)
sorted = sort(in_trial_num);
pre = sorted(diff(sorted) == 1)
post = pre + 1;
end