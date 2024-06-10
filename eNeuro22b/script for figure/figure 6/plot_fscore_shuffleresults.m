load pfc_shuffle_f.mat
f=figure;
h=histogram(shufflef_slope,[-0.04:0.004:0.04]);
h.FaceColor=[94 114 255]/255;
xline(0.0319,'-b','LineWidth',1.5);
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2); %label size
f.Position(3:4)=[600,400]; %figure size
box off;

load ppc_shuffle_f.mat
f=figure;
h=histogram(shufflef_slope,[-0.04:0.004:0.04]);
h.FaceColor=[250 60 60]/255;
xline(0.0142,'-r','LineWidth',1.5);
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2); %label size
f.Position(3:4)=[600,400]; %figure size
box off;