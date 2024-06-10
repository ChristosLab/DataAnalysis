%find significant cells 
[num,txt,raw]=xlsread('neuron_location.xlsx',3);
neuron_num=cell2num([raw(2:end,3)]);
neuron_loc=string([raw(2:end,4)]);
PFC_area={'46','8'};
Parietal_are={'LIP','PPC','7a'};
pfc_index=[];
parietal_index=[];
%locations_list=[15,13,17,5,7]; % 5 locations
location_list=[15,13,11,9,17,1,3,5,7];%from one side to another, 9 locations only match trials
num_physiology_thred=8;
num_behavior_thred=4;
behavior_trials_perblock=40;
trial_touse_edge1=[1,2,3,4,5,6];
trial_touse_edge2=[3,4,5,6,7,8];

for i=1:length(neuron_num)
    if ismember(neuron_loc{i},PFC_area)
    pfc_index=[pfc_index,i];
    else
    parietal_index=[parietal_index,i];
    end
end
pfc_cellnum=neuron_num(pfc_index);
parietal_cellnum=neuron_num(parietal_index);
load('all_msngcorrect_info.mat');
allcell_session=cell_info{1};
load('all_msngcorrect_data.mat');
load('all_behavior_info.mat');
all_session=cell_info{1,1};
for i=1:length(all_session)
    temp_sessionname=all_session{i};
    temp_animalname=temp_sessionname(1:3);
    if strcmpi(temp_animalname,'ken')
        all_animalcode(i)=1;
    else
        all_animalcode(i)=2;
    end
end

behavior_session=behavior_info{1};
behavior_type=behavior_info{2};
behavior_statecode=behavior_info{4};
behavior_reward=behavior_info{5};
%[C,select_index,ib]=intersect(cell_info{3},pfc_cellnum);
[C,select_index,ib]=intersect(cell_info{3},parietal_cellnum);
%select_index=intersect(select_index,find(all_animalcode==2));
prop_category=zeros(length(trial_touse_edge1),4);

for t=1:length(trial_touse_edge1)
good_data=ones(1,length(select_index));
for i=1:length(select_index)
    
    cell_data=msng_data(select_index(i),:); 
    indcell_session=allcell_session{select_index(i)};
    temp_index=find([strcmpi(behavior_session,indcell_session) & (behavior_type==2)]==1);
    if length(temp_index)>0
    cell_statecode=behavior_statecode(temp_index);
    cell_statecode=cell2mat(cell_statecode{1});
    cell_reward=behavior_reward(temp_index);
    cell_reward=cell_reward{1};
    cell_reward=strcmp(cell_reward,'Yes');
    correct_infinish=cell_reward(find(cell_statecode>5));
    
    cell_behavior_blocks(i)=floor(length(correct_infinish)/behavior_trials_perblock);
    cumulative_correct=cumsum(correct_infinish);
    for c=1:cell_behavior_blocks(i)
        correct_inblock=sum(correct_infinish(behavior_trials_perblock*(c-1)+1:behavior_trials_perblock*c));
        correct_proportion(i,c)=correct_inblock/behavior_trials_perblock;
    end
    else    
    correct_proportion(i,:)=NaN(1,length(cell_behavior_blocks(i)));
    end
    
    all_stim=[];
    all_cue=[];
    all_sample=[];
    all_cueloc=[];
    stim_loc=[];
    task_epoch=[];
    for j=1:9 
        cueclass_data=cell_data{location_list(j)}; %use only match trials
        cue_class_trials(j)=length(cueclass_data);
        cueclass_fix_rate=[];
        cueclass_cue_rate=[];
        cueclass_cuedelay_rate=[];
        cueclass_sample_rate=[];
        if cue_class_trials(j)>num_physiology_thred-1
           cueclass_cueon=[cueclass_data.Cue_onT];
           cueclass_sampleon=[cueclass_data.Sample_onT];
           cueclass_TS={cueclass_data.TS};
           for p=1:cue_class_trials(j)
               cueclass_fix_rate(p)=length(find(cueclass_TS{p}>0 & cueclass_TS{p}<cueclass_cueon(p)));
               cueclass_cue_rate(p)=length(find(cueclass_TS{p}>cueclass_cueon(p) & cueclass_TS{p}<cueclass_cueon(p)+0.5))/0.5;
               cueclass_cuedelay_rate(p)=length(find(cueclass_TS{p}>cueclass_cueon(p)+0.5 & cueclass_TS{p}<cueclass_sampleon(p)))/3;
               cueclass_sample_rate(p)=length(find(cueclass_TS{p}>cueclass_sampleon(p) & cueclass_TS{p}<cueclass_sampleon(p)+0.5))/0.5;            
           end
            
   %     if length(cueclass_ismatch)==length(cueclass_fix_count)

        cue_loc=j*ones(1,length(cueclass_cue_rate));
        sample_loc=cue_loc;   
        stim_loc=[stim_loc,cue_loc(trial_touse_edge1(t):trial_touse_edge2(t)),sample_loc(trial_touse_edge1(t):trial_touse_edge2(t))];
        all_cueloc=[all_cueloc,cue_loc(trial_touse_edge1(t):trial_touse_edge2(t))];
        task_epoch=[task_epoch,[ones(1,trial_touse_edge2(t)-trial_touse_edge1(t)+1),2*ones(1,trial_touse_edge2(t)-trial_touse_edge1(t)+1)]];
        all_stim=[all_stim,cueclass_cue_rate(trial_touse_edge1(t):trial_touse_edge2(t)),cueclass_sample_rate(trial_touse_edge1(t):trial_touse_edge2(t))];
        all_cue=[all_cue,cueclass_cue_rate(trial_touse_edge1(t):trial_touse_edge2(t))];
        all_sample=[all_sample,cueclass_sample_rate(trial_touse_edge1(t):trial_touse_edge2(t))];
        temp_cuerate(j)=mean(cueclass_cue_rate(trial_touse_edge1(t):trial_touse_edge2(t)));
        temp_samplerate(j)=mean(cueclass_sample_rate(trial_touse_edge1(t):trial_touse_edge2(t)));
%        end
        end
    end
  
    %location_to_fit=[-2:1:2];
    location_to_fit=[-2,-1,-0.5,-0.25,0,0.25,0.5,1,2];
%{
    try
    gauss_cue(i,:)=coeffvalues(fit(location_to_fit',(temp_cuerate-min(temp_cuerate))','gauss1'));
    gauss_sample(i,:)=coeffvalues(fit(location_to_fit',(temp_samplerate-min(temp_samplerate))','gauss1'));
      if gauss_cue(i,2)>2
         gauss_cue(i,2)=2;
      elseif gauss_cue(i,2)<-2
         gauss_cue(i,2)=-2;
      end
      if gauss_sample(i,2)>2
         gauss_sample(i,2)=2;
      elseif gauss_sample(i,2)<-2
         gauss_sample(i,2)=-2;
      end
    catch
    gauss_cue(i,:)=[0,0,0]; 
    gauss_sample(i,:)=[0,0,0];
    end
    %} 
    if min(cue_class_trials)<num_physiology_thred %filter number of trials
        good_data(i)=0; 
     %   stim_p(i)=[];
     %   delay_p(i)=[];
    end
    cell_numtrials(i)=sum(cue_class_trials);
    
    if length(all_stim)>0 & cue_class_trials(9)>0
    store_cue(i,:)=temp_cuerate;
    store_sample(i,:)=temp_samplerate;
    [temp,raw_cue_preferloc(i)]=max(temp_cuerate);
    [temp,raw_sample_preferloc(i)]=max(temp_samplerate);
    [stim_p(i,:),tbl_stim]=anovan(all_stim,{stim_loc task_epoch},'model','interaction','display','off');
    f_stim(i,:)=cell2mat(tbl_stim(2:4,6))';
    [cue_p(i,:),tbl_cue]=anovan(all_cue,{all_cueloc},'display','off');
    f_cue(i,:)=cell2mat(tbl_cue(2,6));
    [sample_p(i,:),tbl_sample]=anovan(all_sample,{all_cueloc},'display','off');
    f_sample(i,:)=cell2mat(tbl_sample(2,6));
    
    rank_cuerate(i,:)=sort(temp_cuerate/max(temp_cuerate),'descend');
    rank_samplerate(i,:)=sort(temp_samplerate/max(temp_samplerate),'descend');
    store_cuesi(i)=max(temp_cuerate)/((sum(temp_cuerate)-max(temp_cuerate))/(length(location_to_fit)-1));  %selectivity index for each cell
    store_samplesi(i)=max(temp_samplerate)/((sum(temp_samplerate)-max(temp_samplerate))/(length(location_to_fit)-1));
    else
    good_data(i)=0; 
  %  stim_p(i)=[];
  %  delay_p(i)=[];
    end
    
end

    cue_p(cue_p==0) = NaN;
    sample_p(sample_p==0) = NaN;
    stim_p(stim_p==0) = NaN;
    
    data_inclusion_index=find(good_data==1 & cell_behavior_blocks>=num_behavior_thred);
    store_interaction_f(t,:)=f_stim(data_inclusion_index,3);
    behavior_select=correct_proportion(data_inclusion_index,1:num_behavior_thred);

    find_sigstim_main1=intersect(find(stim_p(:,1)<0.05),data_inclusion_index);
    find_sigstim_main2=intersect(find(stim_p(:,2)<0.05),data_inclusion_index);
    find_sigstim_interaction=intersect(find(stim_p(:,3)<0.05),data_inclusion_index);
    behavior_nms_select=correct_proportion(find_sigstim_interaction,1:num_behavior_thred);
    
    find_sigstim_CS1=setdiff(setdiff(find_sigstim_main1,find_sigstim_main2),find_sigstim_interaction);
    find_sigstim_CS2=setdiff(setdiff(find_sigstim_main2,find_sigstim_main1),find_sigstim_interaction);
    find_sigstim_informative=unique([find_sigstim_main1;find_sigstim_main2;find_sigstim_interaction]);
    find_sigstim_LMS=setdiff(intersect(find_sigstim_main1,find_sigstim_main2),find_sigstim_interaction);
    find_sigstim_NMS=find_sigstim_interaction;
    find_sigstim_NMS1=intersect(find_sigstim_interaction,intersect(find_sigstim_main1,find_sigstim_main2));
    find_sigstim_NMS2=intersect(find_sigstim_interaction,setdiff(find_sigstim_main1,find_sigstim_main2));
    find_sigstim_NMS3=intersect(find_sigstim_interaction,setdiff(find_sigstim_main2,find_sigstim_main1));
    find_sigstim_NMS4=setdiff(find_sigstim_interaction,unique([find_sigstim_main1;find_sigstim_main2]));
    find_sigstim_CS=[find_sigstim_CS1;find_sigstim_CS2];
      

    prop_sigstim_main1=length(find_sigstim_main1)/length(data_inclusion_index);
    prop_sigstim_main2=length(find_sigstim_main2)/length(data_inclusion_index);
    prop_sigstim_interaction=length(find_sigstim_interaction)/length(data_inclusion_index);
   %{
    main_effect_y=[prop_sigstim_main1,prop_sigstim_main2,prop_sigstim_interaction];
    bar(main_effect_y);
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
    title('Parietal cue');
    set(gca,'xticklabel',{'stim location','task type','interaction'})
    ylim([0,0.4]);
    text(1,0.4,strcat('n=',num2str(length(data_inclusion_index))));
    %}
    prop_category(t,1)=length(find(cue_p(data_inclusion_index)>=0.05 & sample_p(data_inclusion_index)>=0.05))/length(data_inclusion_index);
    prop_category(t,2)=length(find(cue_p(data_inclusion_index)<0.05 & sample_p(data_inclusion_index)>=0.05))/length(data_inclusion_index);
    prop_category(t,3)=length(find(cue_p(data_inclusion_index)>=0.05 & sample_p(data_inclusion_index)<0.05))/length(data_inclusion_index);
    prop_category(t,4)=length(find(cue_p(data_inclusion_index)<0.05 & sample_p(data_inclusion_index)<0.05))/length(data_inclusion_index);
    meanf_cue(t)=nanmean(f_cue(data_inclusion_index));
    meanf_sample(t)=nanmean(f_sample(data_inclusion_index));
    
    prop_cue_main1(t)=length(find(cue_p(data_inclusion_index)<0.05))/length(data_inclusion_index);
    prop_sample_main1(t)=length(find(sample_p(data_inclusion_index)<0.05))/length(data_inclusion_index);

end
figure;
plot(store_interaction_prop);
title("NMS proportion 2 trial window, 1 trial step" );
xlabel('trial order');
ylabel('Proportion');
c_interaction=polyfit([1:length(store_interaction_prop)],store_interaction_prop*100,1);


figure;
plot(prop_category(:,3));
title("category3 proportion 2 trial window, 1 trial step" );
xlabel('trial order');
ylabel('Proportion');
c_category3=polyfit([1:length(prop_category(:,3))],prop_category(:,3)*100,1);

figure;
plot(prop_cue_main1);
title("cue location selective proportion 2 trial window, 1 trial step" );
xlabel('trial order');
ylabel('Proportion');
c_cueprop=polyfit([1:length(prop_cue_main1)],prop_cue_main1*100,1);

figure;
plot(prop_sample_main1);
title("sample location selective proportion 2 trial window, 1 trial step" );
xlabel('trial order');
ylabel('Proportion');
c_sampleprop=polyfit([1:length(prop_sample_main1)],prop_sample_main1*100,1);

parietal_realresults=c_interaction(1);
save parietal_realresults.mat parietal_realresults
%{
for i=1:size(behavior_select,1)
    behavior_slope(i,:)=polyfit([1:num_behavior_thred],behavior_select(i,1:num_behavior_thred),1);
    interaction_slope(i,:)=polyfit([1:size(store_interaction_f,2)],store_interaction_f(i,:),1);
end
for i=1:size(behavior_nms_select,1)
    behavior_nms_slope(i,:)=polyfit([1:num_behavior_thred],behavior_nms_select(i,1:num_behavior_thred),1);
    interaction_nms_slope(i,:)=polyfit([1:size(store_interaction_f,2)],store_interaction_f(i,:),1);
end
figure;
scatter(behavior_slope(:,1),interaction_slope(:,1));
hold on;
scatter(behavior_nms_slope(:,1),interaction_nms_slope(:,1));
ylim([-2,2]);
title("correlation between behavior and physiological change" );
xlabel('behavior regression slope,4 blocks');
ylabel('NMS proportion regression slope,8 repeat');
legend('all cells','NMS');
%figure;
regression_select=find(isnan(mean(interaction_slope,2))==0);
lm = fitlm(behavior_slope(regression_select,1)',interaction_slope(regression_select,1)');
regression_select=find(isnan(mean(interaction_nms_slope,2))==0);
lm2 = fitlm(behavior_nms_slope(regression_select,1)',interaction_nms_slope(regression_select,1)');
%}
