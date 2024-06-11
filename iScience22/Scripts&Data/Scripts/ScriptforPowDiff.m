function ScriptforPowDiff(Time,Frequency,Spectralpower,TITLE,h)

        T  = Time;  % time events
        f  = Frequency; % frequency
        SS = Spectralpower;% Spectralpower
        pcolor(T,f,SS');shading flat;axis xy;  colormap(h,'polarmap'); colorbar; caxis([-5 5]);

%         contourf(T,f,SS',100,'linecolor','none'); colormap(h,'polarmap');  caxis([-5 5]); colorbar;
        line([0 0],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        line([.5 .5],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        line([2 2],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        line([2.5 2.5],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        line([4 4],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        xlabel("Time (s)"); ylabel("Frequency (Hz)"); axis tight;
        title(TITLE);
        
            
        

