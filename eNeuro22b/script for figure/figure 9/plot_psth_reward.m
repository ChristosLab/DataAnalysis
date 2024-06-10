cell_num_touse=224;
class_touse=[2,3,4,5,6,7];
load('new_msngcorrect_data.mat');
cell_data=msng_data(cell_num_touse,:);


match_class=[15,13,11,9,1,3,5,7];
nonmatch_class=[16,14,12,10,2,4,6,8];
for c=1:length(class_touse)
temp_data=cell_data{match_class(class_touse(c))};
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
plot_end=6;

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
for i=1:length(class_touse)
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
xline(4.5,'LineWidth',1.5);
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2); 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot error trials

load('new_msngerr_data.mat'); %%%%%%%%%%%for err trials
cell_data=msng_data(cell_num_touse,:);
match_class=[15,13,11,9,1,3,5,7];
nonmatch_class=[16,14,12,10,2,4,6,8];
for c=1:length(class_touse)
temp_data=cell_data{nonmatch_class(class_touse(c))};%%%%%%%%%%%for err trials
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
plot_end=6;

color_map(1,:)=[102, 255, 102]/255;
color_map(2,:)=[51, 204, 51]/255;
color_map(3,:)=[0, 153, 51]/255;
color_map(4,:)=[0,51,0]/255;
color_map(5,:)=[0,0,0]/255;
color_map(6,:)=[102,0,102]/255;
color_map(7,:)=[204, 0, 204]/255;
color_map(8,:)=[255, 0, 255]/255;
color_map(9,:)=[255, 102, 255]/255;
subplot(1,2,2);
hold on
for i=1:length(class_touse)
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
ylim([0,temp_ylim(2)]);
xline(4.5,'LineWidth',1.5);
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2); 
