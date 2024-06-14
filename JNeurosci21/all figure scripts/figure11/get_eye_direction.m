%generate fixoff and saccade direction data
%set current path in behavior folder
load('all_spatial_info');
load('all_spatial_data');
mat = dir('*.mat'); 
all_sessions={mat.name};
direct1_classes=[14,16,29,31]; %(10,0)
direct2_classes=[2,4,17,19,34,36];%(10,10)
direct3_classes=[5,7,21,23];%(0,10)
direct4_classes=[9,11,26,28];%(-10,10)
direct5_classes=[13,15,30,32];%(-10,0)
direct6_classes=[1,3,18,20,33,35];%(-10,-10)
direct7_classes=[6,8,22,24];%(0,-10)
direct8_classes=[10,12,25,27];%(10,-10)

for i=1:length(all_sessions)
    [folder, baseFileNameNoExt, extension] = fileparts(all_sessions{i});
    all_sessions_noext{i}=baseFileNameNoExt;
end

for i=1:size(all_spatial_data,1)  %extract task saccade direction from behavior file
    cell_fixoff=[];
    trial_direct=[];
    cell_session=all_spatial_info{i,1};
    if length(cell_session)>8
        cell_session=cell_session(1:8);
    end
    temp_index=find(strcmpi(all_sessions_noext,cell_session));
    if i==2860
        disp('test');
    end
    if length(temp_index)>0
       behavior_data=load(all_sessions{temp_index});
       if isfield(behavior_data.AllData,'trials')
          if length([behavior_data.AllData.trials.FixOff])~=length([behavior_data.AllData.trials.time])
          cell_fixoff=[[behavior_data.AllData.trials.FixOff],6]-[behavior_data.AllData.trials.time];
          else
          cell_fixoff=[behavior_data.AllData.trials.FixOff]-[behavior_data.AllData.trials.time];
          end
          trial_class=[behavior_data.AllData.trials.Class];
          direct1_trials=find(ismember(trial_class,direct1_classes));
          direct2_trials=find(ismember(trial_class,direct2_classes));
          direct3_trials=find(ismember(trial_class,direct3_classes));
          direct4_trials=find(ismember(trial_class,direct4_classes));
          direct5_trials=find(ismember(trial_class,direct5_classes));
          direct6_trials=find(ismember(trial_class,direct6_classes));
          direct7_trials=find(ismember(trial_class,direct7_classes));
          direct8_trials=find(ismember(trial_class,direct8_classes));
          trial_direct(direct1_trials)=1;
          trial_direct(direct2_trials)=2;
          trial_direct(direct3_trials)=3;
          trial_direct(direct4_trials)=4;
          trial_direct(direct5_trials)=5;
          trial_direct(direct6_trials)=6;
          trial_direct(direct7_trials)=7;
          trial_direct(direct8_trials)=8;
          
          all_fixoff{i}=cell_fixoff;
          all_trialdirect{i}=trial_direct;

       end
    end
end
%save spatial_fixoff.mat all_fixoff
save spatial_sacdirect.mat all_trialdirect
load('allcell_post.mat'); 
find_index=allcell_info{10};
%NMS_index=find_index(unique([allcell_info{3};allcell_info{4};allcell_info{5};allcell_info{6}]));
%LS_index=find_index(unique([allcell_info{7};allcell_info{8};allcell_info{9}]));
NMS_index=find_index(unique([allcell_info{3};allcell_info{4};allcell_info{5};allcell_info{6}]));
LS_index=find_index(unique([allcell_info{7};allcell_info{8};allcell_info{9}]));
%find_index=unique([NMS_index;LS_index]); %use all informative
%find_index=LS_index;
find_index=[1736:3040];
cell_count=0;
for i=1:length(find_index)
    cell_data=all_spatial_data(find_index(i),:);
    cell_fixoff=all_fixoff{find_index(i)};
    cell_sacdirect=all_trialdirect{find_index(i)};
    all_cueclass_TS={};
    all_cueclass_trialnum=[];
    all_rewardon=[];
    cue3_TS={};
    cue7_TS={};
    for j=1:8        
       cueclass_data=cell_data{j};
       cueclass_TS={cueclass_data.TS};
       cueclass_cueon={cueclass_data.Cue_onT};
       cueclass_trialnum=[cueclass_data.trialnum];
       cueclass_targon=[cueclass_data.Reward_onT];

       all_rewardon=[all_rewardon,cueclass_targon];
       all_cueclass_TS=[all_cueclass_TS,cueclass_TS];
       if j==3
           cue3_TS=cueclass_TS;
           cue3_cueon=[cueclass_cueon];
       end
       if j==7
           cue7_TS=cueclass_TS;
           cue7_cueon=[cueclass_cueon];
       end
       all_cueclass_trialnum=[all_cueclass_trialnum,cueclass_trialnum];
    end

    for t=1:length(cue3_TS)
        cue3_TS{t}=cue3_TS{t}-cue3_cueon{t};
    end
    for t=1:length(cue7_TS)
        cue7_TS{t}=cue7_TS{t}-cue7_cueon{t};
    end
    
    correct_rewardon=all_rewardon;  %reward time for all correct trials
    allsac_TS={};
    if length(cell_sacdirect)>0
    correct_fixoff=cell_fixoff(all_cueclass_trialnum);
    correct_sacdirect=cell_sacdirect(all_cueclass_trialnum); %sac direction for all correct trials
    cell_count=cell_count+1;
    store_cellindex(cell_count)=find_index(i);
    for p=1:length(all_cueclass_trialnum)
        trial_TS=all_cueclass_TS{p};
        trial_fixoff=correct_fixoff(p);
        allsac_TS{p}=trial_TS-trial_fixoff+0.5;
    end
    direct1_TS=allsac_TS(correct_sacdirect==1);
    direct2_TS=allsac_TS(correct_sacdirect==2);
    all_spikeforchronux=cell2struct(allsac_TS,'spiketime',1);
    direct1_spikeforchronux=cell2struct([direct1_TS],'spiketime',1);
    direct2_spikeforchronux=cell2struct([direct2_TS],'spiketime',1);
    cue3_spikeforchronux=cell2struct([cue3_TS],'spiketime',1);
    cue7_spikeforchronux=cell2struct([cue7_TS],'spiketime',1);
    [psth_all,t_all,E_all,RR_all] = psth(all_spikeforchronux,-0.1,'n');
    [psth_direct1,t_direct1,E_direct1,RR_direct1] = psth(direct1_spikeforchronux,-0.1,'n');
    [psth_direct2,t_direct2,E_direct2,RR_direct2] = psth(direct2_spikeforchronux,-0.1,'n');
    [psth_cue3,t_cue3,E_cue3,RR_cue3] = psth(cue3_spikeforchronux,-0.1,'n');
    [psth_cue7,t_cue7,E_cue7,RR_cue7] = psth(cue7_spikeforchronux,-0.1,'n');
    temp_psth=psth_all(t_all>-3 & t_all<1);
    temp_psth_direct1=psth_direct1(t_direct1>-3 & t_direct1<1);
    temp_psth_direct2=psth_direct2(t_direct2>-3 & t_direct2<1);
    if find_index(i)==1736
        disp('test');
    end
    mean_direct1=max(psth_direct1(t_direct1>-0.5 & t_direct1<0.5));
    mean_direct2=max(psth_direct2(t_direct2>-0.5 & t_direct2<0.5));
    
    maxcue3=max(psth_cue3(t_cue3>-0.2 & t_cue3<0.8));
    maxcue7=max(psth_cue7(t_cue7>-0.2 & t_cue7<0.8));
    %{
    plot(psth_cue3(t_cue3>-1 & t_cue3<5));
    hold on;
    plot(psth_cue7(t_cue3>-1 & t_cue3<5));
    figure;
    plot(psth_all(t_all>-5 & t_all<1));
    hold on;
    plot(psth_direct1(t_direct1>-5 & t_direct1<1));
    plot(psth_direct2(t_direct2>-5 & t_direct2<1));
    close;
%}
    if (length(temp_psth)>=198 & length(temp_psth_direct1)>=198 & length(temp_psth_direct2)>=198)
    store_all_psth(cell_count,:)=temp_psth(1:198);
    norm_all_psth(cell_count,:)=temp_psth(1:198)/max(temp_psth);
    store_direct1_psth(cell_count,:)=temp_psth_direct1(1:198);
    norm_direct1_psth(cell_count,:)=temp_psth_direct1(1:198)/max(temp_psth);
    store_direct2_psth(cell_count,:)=temp_psth_direct2(1:198);
    norm_direct2_psth(cell_count,:)=temp_psth_direct2(1:198)/max(temp_psth);
   %help rhelp r store_all_subcue(cell_count,:)=temp_psth(1:198)-mean([maxcue3,maxcue7]);
   % store_direct1_subcue(cell_count,:)=temp_psth_direct1(1:198)-mean([maxcue3,maxcue7]);
   % store_direct2_subcue(cell_count,:)=temp_psth_direct2(1:198)-mean([maxcue3,maxcue7]);
    tuning_index(cell_count)=abs(mean_direct1-mean_direct2)/(mean_direct1+mean_direct2);
    sac_strength(cell_count)=mean([mean_direct1,mean_direct2]);
    store_cue(cell_count)=mean([maxcue3,maxcue7]);
    store_meanrate(cell_count)=mean(store_all_psth(cell_count,:));
    end
    end
   % store_all_psth{cell_count}=temp_psth(1:74);
    num_timepoint(cell_count)=length(temp_psth);
end
%{
mean_all=mean(store_all_psth);
std_all=1.96*std(store_all_psth)/sqrt(size(store_all_psth,1));
plot_mean_std(mean_all,std_all,'b',0.1,'b',2,'all saccade');
hold on;
mean_direct1=mean(store_direct1_psth);
std_direct1=1.96*std(store_direct1_psth)/sqrt(size(store_direct1_psth,1));
plot_mean_std(mean_direct1,std_direct1,'y',0.1,'y',2,'up saccade');
mean_direct2=mean(store_direct2_psth);
std_direct2=1.96*std(store_direct2_psth)/sqrt(size(store_direct2_psth,1));
plot_mean_std(mean_direct2,std_direct2,'r',0.1,'r',2,'down saccade');
%}
figure;
t=t_all(t_all>-3 & t_all<1);
plot1=plot(t(1:198),mean(store_all_psth),'LineWidth',2);
plot1.Color(4)=0.4;
hold on;
plot2=plot(t(1:198),mean(store_direct1_psth),'LineWidth',2);
plot2.Color(4)=0.4;
plot3=plot(t(1:198),mean(store_direct2_psth),'LineWidth',2);
plot3.Color(4)=0.4;
%hline=refline([0,mean(store_cue)]);
%hline.Color='k';
lgd=legend('all saccade','up saccade','down saccade');
lgd.FontSize=12;
set(lgd,'Position',[0.2,0.65,0.25,0.25]);
%figure;
%plot(mean(norm_all_subcue));
%hold on;
%plot(mean(norm_direct1_subcue));
%plot(mean(norm_direct2_subcue));
%legend('all saccade','up saccade','down saccade');
figure;
scatter(tuning_index,sac_strength);
edges=[0:0.025:0.55];
bincount=histcounts(tuning_index,edges);
figure;
scatter(store_cue(store_meanrate>2),sac_strength(store_meanrate>2));
xlim([0,120]);
ylim([0,120]);
hold on;
plot([0,120],[0,120],'k');
plot([0,60],[0,120],'--k');
set(gca,'fontsize',18)
figure;
edges=[0:16];
hist(sac_strength(store_meanrate>2)./store_cue(store_meanrate>2),edges);
xlim([-1,12]);
%set(gca,'XTickLabel',[]);
%set(gca,'YTickLabel',[]);
set(gca,'fontsize',18)