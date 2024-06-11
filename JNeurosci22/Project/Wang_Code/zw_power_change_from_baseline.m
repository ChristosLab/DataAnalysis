function out_array = zw_power_change_from_baseline(pre_sig, post_sig, f_range, fs)
pre_sp = zw_cmorlet(pre_sig, fs, f_range);
pre_sp = mean(pre_sp, 2);
post_sp = zw_cmorlet(post_sig, fs, f_range);
out_array = (post_sp - pre_sp)./pre_sp;
end