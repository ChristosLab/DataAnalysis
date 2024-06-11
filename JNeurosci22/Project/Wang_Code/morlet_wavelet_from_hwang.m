trialNum = [];
for i = 1:numel(LFPData.class)
    if ~isempty(LFPData.class(i).ntr)
        trialNum = [trialNum, [LFPData.class(i).ntr.Trial_Num]]
    end
end
% tiralNum_as = trialNum
%%
close all
n_cycle = 7;
% Stephanie Jones's wavelet
% g_cmorlet = @(t, f) (f*sqrt(2*pi)/n_cycle).*exp(-(t.^2)*(2*pi^2*f^2)/(n_cycle^2)).*exp(t*(1j*2*pi*f));
% Hwang's wavelet
g_cmorlet = @(t, f) (1/sqrt(2*pi*f)).*exp(-(t.^2)*(2*pi^2*f^2)/(n_cycle^2)).*exp(t*(1j*2*pi*f));

% g_cmorlet = @(lb, ub, n, f) cmorwavf(lb, ub, n, 1/(n_cycle/(2*pi*f))^2, f);
fs = 500;
dur = 2;
rng = 1/fs:1/fs:dur;
target_f_rng = 1:1:60;
f_orig1 = 8;
f_orig2 = 20;
% sig = sin(2*pi*f_orig2*rng);
sig =  sin(2*pi*f_orig1*rng) + cos(2*pi*f_orig2*rng);
% sig = [zeros(1,fs), ones(1, fs), zeros(1,fs)]
% kernel_hw = 10; % Kerel half width in seconds.
kernel_hw = ceil(numel(sig)/fs/2);
kernel_rng = (-kernel_hw):(1/fs):(kernel_hw);
% kernel_lb = -kernel_hw;
% kernel_ub = kernel_hw;
% kernel_n = kernel_hw*fs + 1;
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
x = [rng(1), rng(end)];
y = [target_f_rng(1), target_f_rng(end)];
figure()
image(x, y, out_sig, 'CDataMapping', 'scaled');

% image(x, y, out_sig./max(max(out_sig)), 'CDataMapping', 'scaled');
colorbar();