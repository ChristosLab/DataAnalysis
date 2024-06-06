clear all; close all; clc

% chi square test
n1 = 29; N1 = 151;  %n1 = NS delay; N1 = BS delay
n2 = 63; N2 = 359;  %n2 = NS non-delay; N2 = BS non-delay
x1 = [repmat('a',N1+n1,1); repmat('b',N2+n2,1)]; %x1 = delay column
x2 = [repmat(1,n1,1); repmat(2,N1,1); repmat(1,n2,1); repmat(2,N2,1)]; %x2 = non-delay column
[tbl,chi2stat,pval] = crosstab(x1,x2)

% Number of narrow spiking neurons and total neurons for each group
MSNG_PFC = table([n1;n2],[N1;N2],'VariableNames',{'NS','BS'},'RowNames',{'Delay','Non-delay'})

% Perform Fisher's exact test
[h,p,stats] = fishertest(MSNG_PFC)

%%







