function ScriptforPowS(Time,Frequency,Spectralpower,TITLE,h,numberofTrials)
if nargin < 6; numberofTrials=[]; end;

        T  = Time;  % time events
        f  = Frequency; % frequency
        SS = Spectralpower;% Spectralpower
        pcolor(T,f,SS');shading flat;axis xy;  colormap(h,'parula');  colorbar; caxis([-5 10]);
%         contourf(T,f,SS',100,'linecolor','none'); colormap(h,'parula');  caxis([-5 5]); colorbar;
%         line([-1 -1],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        line([0 0],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        line([.5 .5],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        line([2 2],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        line([2.5 2.5],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        line([4 4],get(h,'Ylim'),'Linestyle',':','LineWidth',1,'Color',[0 0 0]);
        xlabel("Time (s)"); ylabel("Frequency (Hz)"); axis tight;
        if ~isempty(numberofTrials)
            title([TITLE,'-',num2str(numberofTrials)]);
        else
            title(TITLE);
        end

            
        

