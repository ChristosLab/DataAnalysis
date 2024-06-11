function out = zw_ttest2(spectrogram_change_young_ps, spectrogram_change_young_as, spectrogram_change_adult_ps, spectrogram_change_adult_as, ind)
%   Young vs. Adult in PS 
p1 = ttest2(spectrogram_change_young_ps(:, ind), spectrogram_change_adult_ps(:, ind));
%   Young vs. Adult in AS
p2 = ttest2(spectrogram_change_young_as(:, ind), spectrogram_change_adult_as(:, ind));
out = [p1, p2];
end