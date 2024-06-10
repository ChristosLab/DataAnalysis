%cue location vs reward as two factors, include correct match trial and
%incorrect non-match trial
[num,txt,raw]=xlsread('neuron_location.xlsx',3);
neuron_num=cell2num([raw(2:end,3)]);
neuron_loc=string([raw(2:end,4)]);
PFC_area={'46','8'};
Parietal_are={'LIP','PPC','7a'};
pfc_index=[];
parietal_index=[];
location_list_correct=[15,13,11,9,1,3,5,7];%from one side to another, 9 locations
location_list_err=[16,14,12,10,2,4,6,8];
for i=1:length(neuron_num)
    if ismember(neuron_loc{i},PFC_area)
    pfc_index=[pfc_index,i];
    else
    parietal_index=[parietal_index,i];
    end
end
pfc_cellnum=neuron_num(pfc_index);
parietal_cellnum=neuron_num(parietal_index);
load('new_msngerr_info.mat');
load('new_msngcorrect_data.mat');
all_correct_data=msng_data;
load('new_msngerr_data.mat');
all_err_data=msng_data;
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
%animal specific selection

%[C,select_index,ib]=intersect(cell_info{3},pfc_cellnum);
[C,select_index,ib]=intersect(cell_info{3},parietal_cellnum);
%select_index=intersect(select_index,find(all_animalcode==2));
%%%%%%%%%%%%%%%%%%%%store number of trials in correct and error coniditons
for i=1:length(select_index)
cell_data=all_correct_data(select_index(i),:); 
  for j=1:8 
      cueclass_data=cell_data{location_list_correct(j)}; %use only match trials
      store_correct_numtrial(i,j)=length(cueclass_data);
  end
cell_err_data=all_err_data(select_index(i),:); 
  for j=1:8 
      cueclass_data=cell_err_data{location_list_err(j)}; %use only match trials
      store_err_numtrial(i,j)=length(cueclass_data);
  end
cell_err_data=all_err_data(select_index(i),:); 
end
store_correct_numtrial=store_correct_numtrial(:,1:8);%select only match trials
store_err_numtrial=store_err_numtrial(:,1:8);%select only match trials


num_class_threshold=[3,4,5,6];
num_trial_threshold=[2,3,4];
post_filter_count=zeros(length(num_class_threshold),length(num_trial_threshold));
for i=1:length(num_class_threshold)
    %disp(i);
    for j=1:length(num_trial_threshold)
     %   disp(j);
        store_index=[];
        temp_numclass_line=num_class_threshold(i);
        temp_numtrial_line=num_trial_threshold(j);
        for n=1:size(store_err_numtrial,1)
            cell_correct_count=store_correct_numtrial(n,:);
            cell_err_count=store_err_numtrial(n,:);
            if sum((cell_correct_count>=temp_numtrial_line).*(cell_err_count>=temp_numtrial_line))>=temp_numclass_line
            post_filter_count(i,j)=post_filter_count(i,j)+1;
            post_filter_find{i,j,n}=find((cell_correct_count>=temp_numtrial_line).*(cell_err_count>=temp_numtrial_line));
            store_index=[store_index,n];
            end
            post_filter_index{i,j}=store_index;
        end
    end
end

select_class_threshold=1;
select_trial_threshold=3; 
after_filter_index=post_filter_index{select_class_threshold,select_trial_threshold};
after_filter_class=squeeze(post_filter_find(select_class_threshold,select_trial_threshold,after_filter_index));
select_index=select_index(after_filter_index);

numtrial_touse=num_trial_threshold(select_trial_threshold);
good_data=ones(1,length(select_index));
good_err_data=ones(1,length(select_index));

for t=1:50
for i=1:length(select_index)
    cell_class=location_list_correct(after_filter_class{i});
    cell_data=all_correct_data(select_index(i),:); 
    cell_class_err=location_list_err(after_filter_class{i});
    cell_data_err=all_err_data(select_index(i),:);  
    all_reward=[];
    all_cueloc=[];
    if_rewarded=[];
    for j=1:length(cell_class) 
        cueclass_data=cell_data{cell_class(j)}; %use only match trials
        cueclass_data_err=cell_data_err{cell_class_err(j)}; %use only match trials
        cue_class_trials(j)=length(cueclass_data);
        cue_class_trials_err(j)=length(cueclass_data_err);
        trials_touse=randperm(cue_class_trials(j),numtrial_touse);
        trials_touse_err=randperm(cue_class_trials_err(j),numtrial_touse);

        cueclass_reward_rate=[];
        cueclass_reward_rate_err=[];

           cueclass_cueon=[cueclass_data.Cue_onT];
           cueclass_sampleon=[cueclass_data.Sample_onT];
           cueclass_TS={cueclass_data.TS};
           for p=1:cue_class_trials(j)
              cueclass_reward_rate(p)=length(find(cueclass_TS{p}>cueclass_sampleon(p)+1 & cueclass_TS{p}<cueclass_sampleon(p)+2))/1;            
           end
           cueclass_reward_rate=cueclass_reward_rate(trials_touse);

           cueclass_cueon_err=[cueclass_data_err.Cue_onT];
           cueclass_sampleon_err=[cueclass_data_err.Sample_onT];
           cueclass_TS_err={cueclass_data_err.TS};
           for p=1:cue_class_trials_err(j)
              cueclass_reward_rate_err(p)=length(find(cueclass_TS_err{p}>cueclass_sampleon_err(p)+1 & cueclass_TS_err{p}<cueclass_sampleon_err(p)+1.5))/0.5;            
           end
           cueclass_reward_rate_err=cueclass_reward_rate_err(trials_touse_err);
                
         
        cue_loc=j*ones(1,length(cueclass_reward_rate)); 
        all_cueloc=[all_cueloc,cue_loc,cue_loc];
        if_rewarded=[if_rewarded,[ones(1,length(cueclass_reward_rate)),2*ones(1,length(cueclass_reward_rate_err))]];
        all_reward=[all_reward,cueclass_reward_rate,cueclass_reward_rate_err];
    end
   
   if length(all_reward)>0
    [reward_p(i,:),tbl_reward]=anovan(all_reward,{all_cueloc if_rewarded},'model','interaction','display','off');
    f_reward(i,:)=cell2mat(tbl_reward(2:4,6))';
   else
    good_data(i)=0; 
   end
   if isnan(reward_p(i,1))
       good_data(i)=0;
   end
end
    if t==1
    [B,best_reward_index]=sort(reward_p(:,2),'ascend');
    cell_index=select_index(best_reward_index(1:10));
    class_index=after_filter_class(best_reward_index(1:10));
    end
    store_f(:,:,t)=f_reward;

    data_inclusion_index=find(good_data==1);
    find_bothmain=intersect(intersect(find(reward_p(:,1)<0.05),find(reward_p(:,2)<0.05)),data_inclusion_index);   
    find_nomain=intersect(intersect(find(reward_p(:,1)>0.05),find(reward_p(:,2)>0.05)),data_inclusion_index);
    find_sigstim_main1=intersect(find(reward_p(:,1)<0.05),data_inclusion_index);
    find_sigstim_main2=intersect(find(reward_p(:,2)<0.05),data_inclusion_index);
    find_sigstim_interaction=intersect(find(reward_p(:,3)<0.05),data_inclusion_index);
   
  
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
      
    prop_sigstim_main1(t)=length(find_sigstim_main1)/length(data_inclusion_index);
    prop_sigstim_main2(t)=length(find_sigstim_main2)/length(data_inclusion_index);
    prop_sigstim_interaction(t)=length(find_sigstim_interaction)/length(data_inclusion_index);
    
end 
%%%%%%%%%%%%%bar plot for correct and error condition in 50 resamples
figure;
mean_cellf=[nanmean(squeeze(store_f(:,1,:)),2),nanmean(squeeze(store_f(:,2,:)),2),nanmean(squeeze(store_f(:,3,:)),2)];
mean_cellf(find(mean_cellf==Inf))=NaN;
mean_cellf(find(mean_cellf>100))=NaN;
main_effect=nanmean(mean_cellf);
std_resample=std(mean_cellf,0,'omitnan')/sqrt(size(mean_cellf,1));
errhigh_resample=std_resample;
errlow_resample=std_resample;
bar([1,2,3],main_effect);
hold on;
er = errorbar([1,2,3],main_effect,errlow_resample,errhigh_resample);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlimit = get(gca,'xlim');
title('Partietal location vs reward,correct');
set(gca,'xticklabel',{'cue location','if_rewarded','interaction'})


figure;
main_effect_y=[mean(prop_sigstim_main1),mean(prop_sigstim_main2),mean(prop_sigstim_interaction)];
std_y=[std(prop_sigstim_main1)/sqrt(50),std(prop_sigstim_main2)/sqrt(50),std(prop_sigstim_interaction)/sqrt(50)];
errhigh=std_y;
errlow=std_y;
bar([1,2,3],main_effect_y);
hold on;
er = errorbar([1,2,3],main_effect_y,errlow,errhigh);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlimit = get(gca,'xlim');
plot(xlimit,[0.05 0.05]);
title('Partietal cue location vs reward,correct');
set(gca,'xticklabel',{'cue location','if rewarded','interaction'})
%ylim([0,0.3]);
text(1,0.4,strcat('n=',num2str(length(data_inclusion_index))));


