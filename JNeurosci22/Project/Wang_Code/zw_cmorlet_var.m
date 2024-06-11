function out_sig = zw_cmorlet_var(sig, kernel_flag, target_f_range, fs)
%ZW_CMORLET_VAR computes complex Morlet wavelet convolution of input signal
%   kernel_flag - 0: Unscaled wavelet
%               - 1: Jones's wavelet[1]
%               - 2: M.X. Cohen's wavelet[2]
% 
%   target_f_range : Frenquencies to examine in Hertz
%   fs             : Sampling frequency in Hertz

% [1] Jones, S. R., Pritchett, D. L., Sikora, M. A., Stufflebeam, S. M.,
% H�m�l�inen, M., & Moore, C. I. (2009). Quantitative Analysis and
% Biophysically Realistic Neural Modeling of the MEG Mu Rhythm:
% Rhythmogenesis and Modulation of Sensory-Evoked Responses. Journal of
% Neurophysiology, 102(6), 3554�3572. https://doi.org/10.1152/jn.00535.2009
% 
% [2] Cohen, M. X. (2014). Analyzing neural time series data: 
% Theory and practice. The MIT Press.
% 
%%
n_cycle = 7;
if kernel_flag == 0
% Unscaled wavelet
    g_cmorlet = @(t, f) exp(-(t.^2)*(2*pi^2*f^2)/(n_cycle^2)).*exp(t*(1j*2*pi*f));
elseif kernel_flag == 1
    % SJ's wavelet
    g_cmorlet = @(t, f) (f*sqrt(2*pi)/n_cycle).*exp(-(t.^2)*(2*pi^2*f^2)/(n_cycle^2)).*exp(t*(1j*2*pi*f));
elseif kernel_flag == 2
    % M.X. Chohen's wavelet
    g_cmorlet = @(t, f) sqrt(2*f*sqrt(pi)/n_cycle).*exp(-(t.^2)*(2*pi^2*f^2)/(n_cycle^2)).*exp(t*(1j*2*pi*f));
end
% Half width of kernel in seconds
kernel_width = min(1, ceil(numel(sig)/fs/2));
% Points of convolution
kernel_rng = (-kernel_width):(1/fs):(kernel_width);
out_sig = zeros(numel(target_f_range), numel(sig));
ind_f = 0;
for target_f = target_f_range % Loops over frequencies
    ind_f = ind_f + 1;
    %   Actual convolution
    out_sig(ind_f, :) = conv(sig, g_cmorlet(kernel_rng, target_f), 'same');
    %   Convert to density  
    out_sig(ind_f, :) = out_sig(ind_f, :)./fs;
end
%%  Uncomment this part for trial-by-trial visualization
% x = [1, 2];
% y = [target_f_range(1), target_f_range(end)];
% figure()
% out_sig = abs(out_sig).^2;
% image(x, y, out_sig, 'CDataMapping', 'scaled');
% image(x, y, out_sig./max(max(out_sig)), 'CDataMapping', 'scaled');
% colorbar();
end