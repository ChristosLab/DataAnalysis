load('new_msngcorrect_data.mat');
load('parietal_category.mat');
category_save=category_save{4};
%cs1: prefrontal 50 53 54 55 69 81
%cs2: prefrontal 45 51 parietal 60
%LNS: prefrotnal 20  parietal 15 4 5 13
%NMS: prefrontal 8 9 11-14 16-22 29-31 33 parietal 7
%for x=1:length(category_save)
cell_data=msng_data(category_save(7),:);
match_class=[15,13,11,9,17,1,3,5,7];
for c=1:9
temp_data=cell_data{match_class(c)};
cueon=[temp_data.Cue_onT];% cueon time for all trials in the class
spiketimes={temp_data.TS};%spike times for all trials in the class

for i=1:length(cueon)
    temp_TS=spiketimes{i};
    temp_TS=temp_TS-cueon(i);
    match_spiketimes{i}=temp_TS;%cueone aligned spike times for all trials in the struct
end

match_spikeforchronux=cell2struct(match_spiketimes,'spiketime',1);
[psth_match{c},t_match{c},E_match{c},RR] = psth(match_spikeforchronux,0.2,'n'); 
end
plot_start=-1;
plot_end=4;

color_map(1,:)=[102, 255, 102]/255;
color_map(2,:)=[51, 204, 51]/255;
color_map(3,:)=[0, 153, 51]/255;
color_map(4,:)=[0,51,0]/255;
color_map(5,:)=[0,0,0]/255;
color_map(6,:)=[102,0,102]/255;
color_map(7,:)=[204, 0, 204]/255;
color_map(8,:)=[255, 0, 255]/255;
color_map(9,:)=[255, 102, 255]/255;
f_handle=figure;
subplot(1,2,1);
hold on
for i=1:9
    temp_psth=psth_match{i};
    temp_t=t_match{i};
    temp_E=E_match{i};
    cue_rate(i)=mean(temp_psth(temp_t>0  & temp_t<0.5));
    sample_rate(i)=mean(temp_psth(temp_t>3.5  & temp_t<4));
    cue_ratestd(i)=mean(temp_E(temp_t>0  & temp_t<0.5));
    sample_ratestd(i)=mean(temp_E(temp_t>3.5  & temp_t<4));
    curve1=temp_psth+E_match{i};
    curve1=curve1(temp_t>plot_start & temp_t<plot_end);
    curve2=temp_psth-E_match{i};
    curve2=curve2(temp_t>plot_start & temp_t<plot_end);
    inBetween = [curve1, fliplr(curve2)];    
    x_time=linspace(plot_start,plot_end,length(temp_psth(temp_t>plot_start & temp_t<plot_end)));
    x2 = [x_time, fliplr(x_time)];
    fill(x2, inBetween, color_map(i,:),'facealpha',0.1,'LineStyle','none');
    temp_p=plot(x_time,temp_psth(temp_t>plot_start & temp_t<plot_end),'color',color_map(i,:),'LineWidth',1.5);
    cue_count(i)=mean(temp_psth(temp_t>0 & temp_t<0.5));
    cue_err(i)=mean(temp_E(temp_t>0 & temp_t<0.5));
    sample_count(i)=mean(temp_psth(temp_t>3.5 & temp_t<4));
    sample_err(i)=mean(temp_E(temp_t>3.5 & temp_t<4));
    % temp_p.Color(4)=1;

end
temp_ylim=ylim;
ylim([0,temp_ylim(2)]);
%subplot(1,2,2);
f=figure;
p=errorbar([-90,-45,-22.5,-11.25,0,11.25,22.5,45,90],cue_rate,cue_ratestd);
p.LineWidth=1.5;
p.Color=[0,0,204]/255;
hold on;
p=errorbar([-90,-45,-22.5,-11.25,0,11.25,22.5,45,90],sample_rate,sample_ratestd);
p.LineWidth=1.5;
p.Color=[204,0,0]/255;
lgd=legend('cue','sample');
lgd.FontSize = 16;
legend boxoff;
xticks([-90,-45,-22.5,0,22.5,45,90])
xticklabels({'-90','-45','-22.5','0','22.5','45','90'})
xtickangle(-45); %tilt label
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2); %label size
f.Position(3:4)=[600,400]; %figure size
box off;
%xlabel('Stimuli location relative to reference (degrees)');
%ylabel('Firing rate (Hz)');
%ax=gca;
%ax.XAxis.MajorTickChild.LineWidth = 2; % xtick thicker
%ax.YAxis.MajorTickChild.LineWidth = 2; %ytick thicker%
%ax.XAxis.LineWidth=2; %x axis thicker
%ax.YAxis.LineWidth=2;  % yaxis thicker
%end
