function out_sig = zw_cmorlet(sig, fs, target_f_rng)
% tic
rng = 1/fs:1/fs:numel(sig)/fs;
n_cycle = 7;
% Hwang's wavelet
g_cmorlet = @(t, f) (1/sqrt(2*pi*f)).*exp(-(t.^2)*(2*pi^2*f^2)/(n_cycle^2)).*exp(t*(1j*2*pi*f));
% Unscaled wavelet
% g_cmorlet = @(t, f) exp(-(t.^2)*(2*pi^2*f^2)/(n_cycle^2)).*exp(t*(1j*2*pi*f));
kernel_hw = min(1, ceil(numel(sig)/fs/2));
kernel_rng = (-kernel_hw):(1/fs):(kernel_hw);
out_sig = zeros(numel(target_f_rng), numel(sig));
ind_f = 0;
for target_f = target_f_rng
    ind_f = ind_f + 1;
    out_sig(ind_f, :) = abs(...
        conv(sig, g_cmorlet(kernel_rng, target_f), 'same')...
        ).^2;
%     out_sig(ind_f, :) = abs(...
%         conv(sig, g_cmorlet(kernel_lb, kernel_ub, kernel_n, target_f), 'same')...
%         ).^2;
end
% x = [rng(1), rng(end)];
% y = [target_f_rng(1), target_f_rng(end)];
% figure()
% image(x, y, out_sig, 'CDataMapping', 'scaled');
% image(x, y, out_sig./max(max(out_sig)), 'CDataMapping', 'scaled');
% colorbar();
% toc
end