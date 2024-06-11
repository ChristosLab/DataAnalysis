function f = zw_permute_test(set_1, set_2, n_iter)
total_set = [set_1; set_2];
n_1 = numel(set_1);
n_2 = numel(set_2);
resampled = zeros(n_iter, 1);
tic
for i = 1:n_iter
    seq = randperm(n_1 + n_2);
    resampled(i) = mean(total_set(seq(1:n_2))) - mean(total_set(seq(1 + n_2:n_1 + n_2)));
    toc
end
[f, x] = ecdf(resampled);
true_x = mean(set_2) - mean(set_1);
f1 = f(find(true_x > x, 1, 'last'));
f2 = f(find(true_x < x, 1, 'first'));
f = mean([f1, f2]);
if f > 0.5
    f = 1 - f;
end
end