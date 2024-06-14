load('all_spatial_data.mat');
cell_data=all_spatial_data(2884,:);
nonmatch_lookup=[5,6,7,8,1,2,3,4];
for c=1:8
temp_data=cell_data{c};
temp_match_info=[temp_data.IsMatch];
match_index=find(temp_match_info==1);
match_data=temp_data(match_index);
match_cueon=[match_data.Cue_onT];
match_spiketimes={match_data.TS};

temp_data=cell_data{nonmatch_lookup(c)};
temp_match_info=[temp_data.IsMatch];
nonmatch_index=find(temp_match_info==0);
nonmatch_data=temp_data(nonmatch_index);
nonmatch_cueon=[nonmatch_data.Cue_onT];
nonmatch_spiketimes={nonmatch_data.TS};
for i=1:length(match_cueon)
    temp_TS=match_spiketimes{i};
    temp_TS=temp_TS-match_cueon(i);
    match_spiketimes{i}=temp_TS;
end
for i=1:length(nonmatch_cueon)
    temp_TS=nonmatch_spiketimes{i};
    temp_TS=temp_TS-nonmatch_cueon(i);
    nonmatch_spiketimes{i}=temp_TS;
end

match_spikeforchronux=cell2struct(match_spiketimes,'spiketime',1);
nonmatch_spikeforchronux=cell2struct(nonmatch_spiketimes,'spiketime',1);
[psth_match{c},t_match{c},E_match{c},RR] = psth(match_spikeforchronux,0.2,'n'); 
[psth_nonmatch{c},t_nonmatch{c},E_nonmatch{c},RR] = psth(nonmatch_spikeforchronux,0.2,'n'); 
end
%color_map(1,:)=[255,0,0]/255;
%color_map(2,:)=[255,51,153]/255;
%color_map(3,:)=[153,0,255]/255;
%color_map(4,:)=[0,0,255]/255;
%color_map(5,:)=[0,204,255]/255;
%color_map(6,:)=[0,255,0]/255;
%color_map(7,:)=[255,255,0]/255;
%color_map(8,:)=[255,102,0]/255;

color_map(1,:)=[243,233,28]/255;
color_map(2,:)=[147,219,53]/255;
color_map(3,:)=[61,195,108]/255;
color_map(4,:)=[34,162,135]/255;
color_map(5,:)=[52,127,142]/255;
color_map(6,:)=[67,92,141]/255;
color_map(7,:)=[79,48,127]/255;
color_map(8,:)=[72,0,84]/255;
f_handle=figure;
hold on
for i=1:8
    temp_psth=psth_match{i};
    temp_t=t_match{i};
    temp_E=E_match{i};
    curve1=temp_psth+E_match{i};
    curve1=curve1(temp_t>0 & temp_t<4);
    curve2=temp_psth-E_match{i};
    curve2=curve2(temp_t>0 & temp_t<4);
    inBetween = [curve1, fliplr(curve2)];    
    x_time=linspace(0,4,length(temp_psth(temp_t>0 & temp_t<4)));
    x2 = [x_time, fliplr(x_time)];
    fill(x2, inBetween, color_map(i,:),'facealpha',0.1,'LineStyle','none');
    temp_p=plot(x_time,temp_psth(temp_t>0 & temp_t<4),'color',color_map(i,:),'LineWidth',1.5);
    match_count(i)=mean(temp_psth(temp_t>2 & temp_t<2.5));
    match_err(i)=mean(temp_E(temp_t>2 & temp_t<2.5));
    % temp_p.Color(4)=1;
    temp_psth=psth_nonmatch{i};
    temp_t=t_nonmatch{i};
    temp_E=E_nonmatch{i};
    curve1=temp_psth+E_nonmatch{i};
    curve1=curve1(temp_t>0 & temp_t<4);
    curve2=temp_psth-E_nonmatch{i};
    curve2=curve2(temp_t>0 & temp_t<4);
    inBetween = [curve1, fliplr(curve2)];  
    x_time=linspace(0,4,length(temp_psth(temp_t>0 & temp_t<4)));
    x2 = [x_time, fliplr(x_time)];
    fill(x2, inBetween, color_map(i,:),'facealpha',0.1,'LineStyle','none');
    temp_p=plot(x_time,temp_psth(temp_t>0 & temp_t<4),'--','color',color_map(i,:),'LineWidth',1.5);
    nonmatch_count(i)=mean(temp_psth(temp_t>2 & temp_t<2.5));
    nonmatch_err(i)=mean(temp_E(temp_t>2 & temp_t<2.5));
    % temp_p.Color(4)=0.15;
end
temp_ylim=ylim;
ylim([0,temp_ylim(2)]);
fill([2,2.5,2.5,2],[0,0,temp_ylim(2),temp_ylim(2)],[128,128,128]/255,'facealpha',0.3,'LineStyle','none');


