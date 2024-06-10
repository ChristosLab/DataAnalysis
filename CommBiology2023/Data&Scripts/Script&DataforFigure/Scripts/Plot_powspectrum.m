function Plot_powspectrum(Time,Freqs,Pow,TIT,h)
pcolor(Time,Freqs,10*log10(Pow)');shading flat;axis xy;  colormap(h,'parula'); caxis([0 2]); colorbar;
line([0 0],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
line([.5 .5],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
line([2 2],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
line([2.5 2.5],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
line([4 4],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
xlabel("Time (s)"); ylabel("Frequency (Hz)"); axis tight; ylim([0 128]);
if ~isempty(TIT)
    title(TIT);
end
    

            
        
