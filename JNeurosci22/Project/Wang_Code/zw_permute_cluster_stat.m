function [bw, corrected_p] = zw_permute_cluster_stat(stat, p, critical_p, temp_input, n_permute, stat_table_order)
[bw, n] = bwlabel(p < critical_p);
input = cat(1, temp_input{:});
size_ = [size(temp_input{1}, 1), size(temp_input{2}, 1), size(temp_input{3}, 1), size(temp_input{4}, 1)];
range_ = [cumsum([1, size_(1:end-1)])', cumsum(size_)'];
corrected_p = zeros(n, 1);
for i = 1:n
    tic
    elements_ = (bw == i);
    input_ = input(:, elements_);
    orig_stat_ = sum(stat(elements_));
    perm_stat_ = zeros(n_permute, 1);
    parfor j = 1:n_permute
        permed_ind = randperm(size(input, 1));
        input__ = {...
            input_(permed_ind(range_(1,1):range_(1,2)), :,:),...
            input_(permed_ind(range_(2,1):range_(2,2)), :,:),...
            input_(permed_ind(range_(3,1):range_(3,2)), :,:),...
            input_(permed_ind(range_(4,1):range_(4,2)), :,:)...
            };
        tbl_ = zw_group_anova(input__);
        perm_stat_(j) = tbl_{stat_table_order, 6, :};
    end
    [cdf_f, cdf_x] = ecdf(perm_stat_);
    corrected_p(i) = mean(...
        cdf_f([find(orig_stat_ > cdf_x, 1), find(orig_stat_ < cdf_x, 1)])...
        );
    toc
    sprintf('One cluster done')
end
end