stft_y = [200,1100];
cwt_y = [1000, 16000];
%%
subj_to_plot = 2;
%%
figure();
imagesc(t_cwt, [2, 100],mean(cwt_young_ps{subj_to_plot}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
figure();
imagesc(t_cwt, [2, 100],mean(cwt_adult_ps{subj_to_plot}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
%%
figure();
imagesc(t_stft, [2, 100],mean(log(stft_young_ps{subj_to_plot}(2:end, :, :)), 3), log(stft_y));set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
figure();
imagesc(t_stft, [2, 100],mean(log(stft_adult_ps{subj_to_plot}(2:end, :, :)), 3), log(stft_y));set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
%%
figure();
imagesc(t_stft, [2, 100],mean(stft_young_ps{subj_to_plot}(2:end, :, :), 3), stft_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
figure();
imagesc(t_stft, [2, 100],mean(stft_adult_ps{subj_to_plot}(2:end, :, :), 3), stft_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
%%
% 
% 
% 
% 
%%
figure();
imagesc(t_cwt, [2, 100],mean(cwt_young_as{subj_to_plot}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
figure();
imagesc(t_cwt, [2, 100],mean(cwt_adult_as{subj_to_plot}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
%%
figure();
imagesc(t_stft, [2, 100],mean(log(stft_young_as{subj_to_plot}(2:end, :, :)), 3), log(stft_y));set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
figure();
imagesc(t_stft, [2, 100],mean(log(stft_adult_as{subj_to_plot}(2:end, :, :)), 3), log(stft_y));set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
%%
figure();
imagesc(t_stft, [2, 100],mean(stft_young_as{subj_to_plot}(2:end, :, :), 3), stft_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
figure();
imagesc(t_stft, [2, 100],mean(stft_adult_as{subj_to_plot}(2:end, :, :), 3), stft_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()