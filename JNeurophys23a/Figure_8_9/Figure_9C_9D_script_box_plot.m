BrainArea='PPC';
No_neurons=147;
figure
% subplot 221
% 
% tbl=readtable('distance_table.xlsx','Sheet',[BrainArea, '_cue']);
% Order = {'R1','R2'};
% tbl.R1R2 = categorical(tbl.R1R2,Order);
% boxchart(tbl.R1R2,tbl.Distance,'GroupByColor',tbl.CorrErr)
% ylabel('Distance')
% legend
% title([BrainArea, ' cue'])
% 
% clear tbl Order

% subplot 222

tbl=readtable('distance_table.xlsx','Sheet',[BrainArea, '_cd']);
Order = {'Remember 1st','Remember 2nd'};
tbl.R1R2 = categorical(tbl.R1R2,Order);
boxchart(tbl.R1R2,tbl.Distance,'GroupByColor',tbl.CorrErr)
ylabel('Distance')
legend
title(sprintf('PPC cuedelay, n=%d', No_neurons))

clear tbl Order

% subplot 223
% 
% tbl=readtable('distance_table.xlsx','Sheet', [BrainArea, '_sample']);
% Order = {'R1','R2'};
% tbl.R1R2 = categorical(tbl.R1R2,Order);
% boxchart(tbl.R1R2,tbl.Distance,'GroupByColor',tbl.CorrErr)
% ylabel('Distance')
% legend
% title([BrainArea, ' sample'])
% 
% clear tbl Order
% 
% 
% subplot 224
% 
% tbl=readtable('distance_table.xlsx', 'Sheet',[BrainArea, '_sd']);
% Order = {'R1','R2'};
% tbl.R1R2 = categorical(tbl.R1R2,Order);
% boxchart(tbl.R1R2,tbl.Distance,'GroupByColor',tbl.CorrErr)
% ylabel('Distance')
% legend
% title([BrainArea, ' sample delay            '])
% 
% clear tbl Order
% 
% sgtitle('PPC Cue')

% sgtitle([BrainArea])






