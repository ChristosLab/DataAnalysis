function psth = zw_spike_time_to_psth(st, bin_width, step, t_range)
fine_step_scale = gcd(1000*bin_width, 1000*step);
fine_step = fine_step_scale/1000;
ratio_bin = 1000*bin_width/fine_step_scale;
ratio_step = 1000*step/fine_step_scale;
bin_edges = t_range(1):fine_step:t_range(2);
windows = (0:ratio_bin - 1) + (1:ratio_step:(numel(bin_edges) - 1))';
bin_edges = (t_range(1) - (bin_width - step)/2):fine_step:(t_range(2) + (bin_width - step)/2);
fine_spike_count = histcounts(st, bin_edges);
psth = sum(fine_spike_count(windows), 2)'/bin_width;
end