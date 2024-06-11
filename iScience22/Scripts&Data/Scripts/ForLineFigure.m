function ForLineFigure(Data,Time,CL)
D          = Data;
T          = Time;
data_mean  = mean(D,1); % mean &error bars
SORTPS     = sort(D,1);
lE         = SORTPS(50,:);
uE         = SORTPS(950,:);
yP         = [lE,fliplr(uE)];%Make the patch & remove nans 
xP         = [T,fliplr(T)];
yP(isnan(yP)) = [];

lineProps  = {CL, 'LineWidth', 1};
H.mainLine = plot(T,data_mean,lineProps{:}); hold on
col        = get(H.mainLine,'color');
edgeColor  = col+(1-col)*0.6;
faceAlpha  = 0.05;                              
patchColor = col;
H.patch    = patch(xP,yP,1,'facecolor',patchColor,'edgecolor','none','facealpha',faceAlpha);
H.edge(1)  = plot(T,lE,'-','color',edgeColor);
H.edge(2)  = plot(T,uE,'-','color',edgeColor);
H.mainLine = plot(T,data_mean,lineProps{:});
clear yP xP lE uE SORTPS data_mean D Time CL;
