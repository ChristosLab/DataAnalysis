function psth_plot(Neuronaldata, t,Smoothing,linetypes,colormap,Label)
yl = [1, 30];

hold on;
for jj=1:size(Neuronaldata,2)
    PSTH  = Neuronaldata{jj};
    if strcmpi(Smoothing,'y')
        windowWidth = 10;
        averagedData = zeros(size(PSTH,1), size(PSTH,2) - windowWidth);
        for k = 1 : size(PSTH,2) - windowWidth
            averagedData(:,k) = mean(PSTH(:,k:k+windowWidth-1),2);
        end
        t1=t(1):0.02:t(end-windowWidth);
        if size(averagedData,1)==1
            p(jj)=plot(t1,averagedData,'linestyle',linetypes{jj},'color',colormap{jj},'LineWidth', 1);
        else
            fill([t1, fliplr(t1)],[nanmean(averagedData, 1) + nanstd(averagedData, 1)/sqrt(size(averagedData,1)),...
                fliplr(nanmean(averagedData, 1) - nanstd(averagedData, 1)/sqrt(size(averagedData,1)))], colormap{jj}, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
            p(jj)=plot(t1,nanmean(averagedData),'linestyle',linetypes{jj},'color',colormap{jj},'LineWidth',1);
            
        end
    else
        if size(PSTH,1)==1
            p(jj)=plot(t,PSTH,'linestyle',linetypes{jj},'color',colormap{jj},'LineWidth', 1);
        else
           fill([t, fliplr(t)],[nanmean(PSTH, 1) + nanstd(PSTH, 1)/sqrt(size(PSTH,1)),...
                fliplr(nanmean(PSTH, 1) - nanstd(PSTH, 1)/sqrt(size(PSTH,1)))], colormap{jj}, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
            p(jj)=plot(t,nanmean(PSTH),'linestyle',linetypes{jj},'color',colormap{jj},'LineWidth',1);
        end
    end
  
end


if isempty(yl)
    yl = ylim;
end
ylim(yl);
line([0 0],yl,'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
line([.5 .5],yl,'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
line([2 2],yl,'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
line([2.5 2.5],yl,'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
line([4 4],yl,'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
ylabel('Firing rate (sp/s)'); xlabel("Time (s)"); 

if strcmpi(Smoothing,'y')
    axis tight; xlim([t1(1) t1(end)]); 
else
    axis tight; xlim([t(1) t(end)]);
end

 legend(p,{Label{:}}); 
