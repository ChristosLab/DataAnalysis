function out = zw_pev(anova_tbl)
%ZW_PEV computes the unbiased porpotion of variance explained from a MATLAB
%anova output
df_group = anova_tbl{2,3};
ss_total = anova_tbl{4, 2};
ss_group = anova_tbl{2, 2};
ms_error = anova_tbl{3, 4};
out = (ss_group - df_group*ms_error)/(ss_total + ms_error);
end