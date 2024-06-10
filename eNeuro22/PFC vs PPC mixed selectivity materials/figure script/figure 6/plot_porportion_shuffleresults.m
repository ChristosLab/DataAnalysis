load pfc_shuffleresults.mat
f=figure;
h=histogram(shuffle_results(1,:),[-0.6:0.05:0.6]);
h.FaceColor=[94 114 255]/255;
load pfc_realresults.mat
xline(pfc_realresults,'-b','LineWidth',1.5);
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2); %label size
f.Position(3:4)=[600,400]; %figure size
box off;

load parietal_shuffleresults.mat
f=figure;
h=histogram(shuffle_results(1,:),[-0.6:0.05:0.6]);
h.FaceColor=[250 60 60]/255;
load parietal_realresults.mat
xline(parietal_realresults,'-r','LineWidth',1.5);
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2); %label size
f.Position(3:4)=[600,400]; %figure size
box off;