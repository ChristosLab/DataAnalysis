function out = zw_pac_extract(repo, ind, low_freq, hi_freq)
out = zeros(numel(ind), 1);
for i = 1:numel(ind)
    temp_out = [];
    for j = 1:8
    temp_out = [temp_out, mean(repo(ind(i)).class(j).normalized_pac(1, hi_freq, low_freq), 'all')];
    end
    out(i) = mean(temp_out);
end
end