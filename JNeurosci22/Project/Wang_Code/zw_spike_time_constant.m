function [tau, exit_code, rho_mean, ci] = zw_spike_time_constant(mat, bin_width, max_lag, average_time_points)
%ZW_SPIKE_TIME_CONSTANT estimates the innate decay time constant from a
%matrix of spike trains
%   Inputs:
%          mat - N trial X Time point firing rate
%          bin_width - milliscond length of nonoverlapping bins
%          Maximum distance of time point pairs for correlation
firng_rate_thresh = 1;
first_rho_decrease_time = 150;
A_bounds   = [0, inf];
tau_bounds = [0, 500];
first_rho_decrease_point = round(first_rho_decrease_time/bin_width);

len = size(mat, 2);
n_delta = max_lag/bin_width;
delta = (1:n_delta)*bin_width; %    Lags in milliseconds
rho_mat = corrcoef(mat);
rho      = cell(1, n_delta);
rho_mean = nan(1, n_delta);
ci = nan(3, 2);
if mean(mat, 'all') < firng_rate_thresh
    tau = nan;
    exit_code = 0;
    disp('Mean firing rate too low')
    return
end
for i = 1:n_delta
    rho_ = zeros(1, len - i);
    for j = 1:(len - i)
        rho_(j) = rho_mat(j, j + i);
    end
    rho{i} = rho_;
    rho_mean(i) = nanmean(rho_);
end
% if any(diff(rho(1:(first_rho_decrease_point - 1))) < 0)
%     tau = nan;
%     exit_code = -1;
%     disp('Rho decreasing at too small lag')
%     rho
%     return
% end
first_decrease = find(diff(rho_mean) < 0, 1);
if average_time_points
    delta_fin = delta(first_decrease:end);
    rho_fin   = rho_mean(first_decrease:end);
else
    delta_ = [];
    for i = first_decrease:n_delta
    delta_ = [delta_, zeros(size(rho{i})) + delta(i)];
    end
    delta_fin = delta_;
    rho_fin   = [rho{first_decrease:end}];
end
%   para: 1 - A; 2 - tau; 3 - B
%   rho(delta) = A(exp(-delta/tau) + B)
min_sq = @(para) rho_fin - para(1).*(exp(-delta_fin/para(2)) + para(3));
options.Algorithm = 'levenberg-marquardt';
% options.Display   = 'off';
[out_para, ~, res, ~, ~, ~, jacobian]  = lsqnonlin(min_sq, [1, 300, 1], [], [], options);
if out_para(1) < A_bounds(1) || out_para(1) > A_bounds(2)
    tau       = nan;
    exit_code = -2;
    fprintf('A = %.2f out of range\n', out_para(1));
    return
elseif out_para(2) < tau_bounds(1) || out_para(2) > tau_bounds(2)
    tau       = nan;
    exit_code = -3;
    disp('Tau out of range')
    return
end
tau = out_para(2);
exit_code = 1;
ci = nlparci(out_para, res, 'jacobian', jacobian);
end