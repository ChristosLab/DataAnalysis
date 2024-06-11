function out = zw_pev_v2(anova_tbl)
%ZW_PEV computes the unbiased porpotion of variance explained from a MATLAB
%anova output V2: Updated for cases of zero squaured error (arising from
%equal firing across all trials at a time point) PEV is set to 0 in this
%situation 
%2/7/21 - ZW
df_group = anova_tbl{2,3};
ss_total = anova_tbl{4, 2};
ss_group = anova_tbl{2, 2};
ms_error = anova_tbl{3, 4};
out = (ss_group - df_group*ms_error)/(ss_total + ms_error);
if isnan(out)
    out = 0;
end
end