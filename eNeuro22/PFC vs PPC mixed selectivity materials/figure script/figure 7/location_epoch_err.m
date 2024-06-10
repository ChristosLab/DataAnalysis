%find significant cells 
[num,txt,raw]=xlsread('neuron_location.xlsx',3);
neuron_num=cell2num([raw(2:end,3)]);
neuron_loc=string([raw(2:end,4)]);
PFC_area={'46','8'};
Parietal_are={'LIP','PPC','7a'};
pfc_index=[];
parietal_index=[];
%locations_list=[15,13,17,5,7]; % 5 locations
location_list=[15,13,11,9,17,1,3,5,7];%from one side to another, 9 locations
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

[C,select_index,ib]=intersect(cell_info{3},pfc_cellnum);%
%[C,select_index,ib]=intersect(cell_info{3},parietal_cellnum);
%select_index=intersect(select_index,find(all_animalcode==2));
%%%%%%%%%%%%%%%%%%%%store number of trials in correct and error coniditons
for i=1:length(select_index)
cell_data=all_correct_data(select_index(i),:); 
  for j=1:9 
      cueclass_data=cell_data{location_list(j)}; %use only match trials
      store_correct_numtrial(i,j)=length(cueclass_data);
  end
cell_err_data=all_err_data(select_index(i),:); 
  for j=1:9 
      cueclass_data=cell_err_data{location_list(j)}; %use only match trials
      store_err_numtrial(i,j)=length(cueclass_data);
  end
cell_err_data=all_err_data(select_index(i),:); 
end
store_correct_numtrial=store_correct_numtrial(:,1:9);%select only match trials
store_err_numtrial=store_err_numtrial(:,1:9);%select only match trials


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
    cell_class=location_list(after_filter_class{i});
    cell_data=all_correct_data(select_index(i),:);  
    all_stim=[];
    all_cue=[];
    all_sample=[];
    all_cueloc=[];
    stim_loc=[];
    task_epoch=[];
    for j=1:length(cell_class) 
        cueclass_data=cell_data{cell_class(j)}; %use only match trials
        cue_class_trials(j)=length(cueclass_data);
        trials_touse=randperm(cue_class_trials(j),numtrial_touse);
        cueclass_fix_rate=[];
        cueclass_cue_rate=[];
        cueclass_cuedelay_rate=[];
        cueclass_sample_rate=[];

           cueclass_cueon=[cueclass_data.Cue_onT];
           cueclass_sampleon=[cueclass_data.Sample_onT];
           cueclass_TS={cueclass_data.TS};
           for p=1:cue_class_trials(j)
               cueclass_fix_rate(p)=length(find(cueclass_TS{p}>0 & cueclass_TS{p}<cueclass_cueon(p)));
               cueclass_cue_rate(p)=length(find(cueclass_TS{p}>cueclass_cueon(p) & cueclass_TS{p}<cueclass_cueon(p)+0.5))/0.5;
               cueclass_cuedelay_rate(p)=length(find(cueclass_TS{p}>cueclass_cueon(p)+0.5 & cueclass_TS{p}<cueclass_sampleon(p)))/3;
               cueclass_sample_rate(p)=length(find(cueclass_TS{p}>cueclass_sampleon(p) & cueclass_TS{p}<cueclass_sampleon(p)+0.5))/0.5;            
           end
         cueclass_fix_rate=cueclass_fix_rate(trials_touse); 
         cueclass_cue_rate=cueclass_cue_rate(trials_touse);
         cueclass_cuedelay_rate=cueclass_cuedelay_rate(trials_touse);
         cueclass_sample_rate=cueclass_sample_rate(trials_touse);
   %     if length(cueclass_ismatch)==length(cueclass_fix_count)

        cue_loc=j*ones(1,length(cueclass_cue_rate));
        sample_loc=cue_loc;   
        stim_loc=[stim_loc,cue_loc,sample_loc];
        all_cueloc=[all_cueloc,cue_loc];
        task_epoch=[task_epoch,[ones(1,length(cueclass_cue_rate)),2*ones(1,length(cueclass_cue_rate))]];
        all_stim=[all_stim,cueclass_cue_rate,cueclass_sample_rate];
        all_cue=[all_cue,cueclass_cue_rate];
        all_sample=[all_sample,cueclass_sample_rate];

    end

    
   if length(all_stim)>0
    [stim_p(i,:),tbl_stim]=anovan(all_stim,{stim_loc task_epoch},'model','interaction','display','off');
    f_stim(i,:)=cell2mat(tbl_stim(2:4,6))';
    [cue_p(i,:),tbl_cue]=anovan(all_cue,{all_cueloc},'display','off');
    f_cue(i,:)=cell2mat(tbl_cue(2,6));
    [sample_p(i,:),tbl_sample]=anovan(all_sample,{all_cueloc},'display','off');
    f_sample(i,:)=cell2mat(tbl_sample(2,6)); 
   else
    good_data(i)=0; 
   end

end



for i=1:length(select_index)
    cell_class=location_list(after_filter_class{i});
    cell_data=all_err_data(select_index(i),:);  
    all_stim=[];
    all_cue=[];
    all_sample=[];
    all_cueloc=[];
    stim_loc=[];
    task_epoch=[];
    for j=1:length(cell_class) 
        cueclass_data=cell_data{cell_class(j)}; %use only match trials
        cue_class_trials(j)=length(cueclass_data);
        trials_touse=randperm(cue_class_trials(j),numtrial_touse);
        cueclass_fix_rate=[];
        cueclass_cue_rate=[];
        cueclass_cuedelay_rate=[];
        cueclass_sample_rate=[];

           cueclass_cueon=[cueclass_data.Cue_onT];
           cueclass_sampleon=[cueclass_data.Sample_onT];
           cueclass_TS={cueclass_data.TS};
           for p=1:cue_class_trials(j)
               cueclass_fix_rate(p)=length(find(cueclass_TS{p}>0 & cueclass_TS{p}<cueclass_cueon(p)));
               cueclass_cue_rate(p)=length(find(cueclass_TS{p}>cueclass_cueon(p) & cueclass_TS{p}<cueclass_cueon(p)+0.5))/0.5;
               cueclass_cuedelay_rate(p)=length(find(cueclass_TS{p}>cueclass_cueon(p)+0.5 & cueclass_TS{p}<cueclass_sampleon(p)))/3;
               cueclass_sample_rate(p)=length(find(cueclass_TS{p}>cueclass_sampleon(p) & cueclass_TS{p}<cueclass_sampleon(p)+0.5))/0.5;            
           end
         cueclass_fix_rate=cueclass_fix_rate(trials_touse); 
         cueclass_cue_rate=cueclass_cue_rate(trials_touse);
         cueclass_cuedelay_rate=cueclass_cuedelay_rate(trials_touse);
         cueclass_sample_rate=cueclass_sample_rate(trials_touse);
   %     if length(cueclass_ismatch)==length(cueclass_fix_count)

        cue_loc=j*ones(1,length(cueclass_cue_rate));
        sample_loc=cue_loc;   
        stim_loc=[stim_loc,cue_loc,sample_loc];
        all_cueloc=[all_cueloc,cue_loc];
        task_epoch=[task_epoch,[ones(1,length(cueclass_cue_rate)),2*ones(1,length(cueclass_cue_rate))]];
        all_stim=[all_stim,cueclass_cue_rate,cueclass_sample_rate];
        all_cue=[all_cue,cueclass_cue_rate];
        all_sample=[all_sample,cueclass_sample_rate];

    end

    
   if length(all_stim)>0
    [stim_err_p(i,:),tbl_stim]=anovan(all_stim,{stim_loc task_epoch},'model','interaction','display','off');
    f_err_stim(i,:)=cell2mat(tbl_stim(2:4,6))';
    [cue_err_p(i,:),tbl_cue]=anovan(all_cue,{all_cueloc},'display','off');
    f_err_cue(i,:)=cell2mat(tbl_cue(2,6));
    [sample_err_p(i,:),tbl_sample]=anovan(all_sample,{all_cueloc},'display','off');
    f_err_sample(i,:)=cell2mat(tbl_sample(2,6)); 
   else
    good_err_data(i)=0; 
   end
end
    store_f_correct(:,:,t)=f_stim;
    store_f_err(:,:,t)=f_err_stim;
    
    cue_p(cue_p==0) = NaN;
    sample_p(sample_p==0) = NaN;
    stim_p(stim_p==0) = NaN;
    
    data_inclusion_index=find(good_data==1);
    find_bothmain=intersect(intersect(find(stim_p(:,1)<0.05),find(stim_p(:,2)<0.05)),data_inclusion_index);   
    find_nomain=intersect(intersect(find(stim_p(:,1)>0.05),find(stim_p(:,2)>0.05)),data_inclusion_index);
    find_sigstim_main1=intersect(find(stim_p(:,1)<0.05),data_inclusion_index);
    find_sigstim_main2=intersect(find(stim_p(:,2)<0.05),data_inclusion_index);
    find_sigstim_interaction=intersect(find(stim_p(:,3)<0.05),data_inclusion_index);
   
  
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
    
   
   
   %................err proportion................. 
    cue_err_p(cue_err_p==0) = NaN;
    sample_err_p(sample_err_p==0) = NaN;
    stim_err_p(stim_err_p==0) = NaN;
    
    data_inclusion_err_index=find(good_err_data==1);
    find_bothmain_err=intersect(intersect(find(stim_err_p(:,1)<0.05),find(stim_err_p(:,2)<0.05)),data_inclusion_err_index);   
    find_nomain_err=intersect(intersect(find(stim_err_p(:,1)>0.05),find(stim_err_p(:,2)>0.05)),data_inclusion_err_index);
    find_sigstim_main1_err=intersect(find(stim_err_p(:,1)<0.05),data_inclusion_err_index);
    find_sigstim_main2_err=intersect(find(stim_err_p(:,2)<0.05),data_inclusion_err_index);
    find_sigstim_interaction_err=intersect(find(stim_err_p(:,3)<0.05),data_inclusion_err_index);
   
  
    find_sigstim_CS1_err=setdiff(setdiff(find_sigstim_main1_err,find_sigstim_main2_err),find_sigstim_interaction_err);
    find_sigstim_CS2_err=setdiff(setdiff(find_sigstim_main2_err,find_sigstim_main1_err),find_sigstim_interaction_err);
    find_sigstim_informative_err=unique([find_sigstim_main1_err;find_sigstim_main2_err;find_sigstim_interaction_err]);
    find_sigstim_LMS_err=setdiff(intersect(find_sigstim_main1_err,find_sigstim_main2_err),find_sigstim_interaction_err);
    find_sigstim_NMS_err=find_sigstim_interaction_err;
    find_sigstim_NMS1_err=intersect(find_sigstim_interaction_err,intersect(find_sigstim_main1_err,find_sigstim_main2_err));
    find_sigstim_NMS2_err=intersect(find_sigstim_interaction_err,setdiff(find_sigstim_main1_err,find_sigstim_main2_err));
    find_sigstim_NMS3_err=intersect(find_sigstim_interaction_err,setdiff(find_sigstim_main2_err,find_sigstim_main1_err));
    find_sigstim_NMS4_err=setdiff(find_sigstim_interaction_err,unique([find_sigstim_main1_err;find_sigstim_main2_err]));
    find_sigstim_CS_err=[find_sigstim_CS1;find_sigstim_CS2];
      

    prop_sigstim_err_main1(t)=length(find_sigstim_main1_err)/length(data_inclusion_err_index);
    prop_sigstim_err_main2(t)=length(find_sigstim_main2_err)/length(data_inclusion_err_index);
    prop_sigstim_err_interaction(t)=length(find_sigstim_interaction_err)/length(data_inclusion_err_index);
    
end 

figure;
mean_cellf_correct=[nanmean(squeeze(store_f_correct(:,1,:)),2),nanmean(squeeze(store_f_correct(:,2,:)),2),nanmean(squeeze(store_f_correct(:,3,:)),2)];
mean_cellf_correct(find(mean_cellf_correct==Inf))=NaN;
mean_cellf_correct(find(mean_cellf_correct>100))=NaN;
main_effect_correct=nanmean(mean_cellf_correct);
std_correct=std(mean_cellf_correct,0,'omitnan')/sqrt(size(mean_cellf_correct,1));
errhigh_correct=std_correct;
errlow_correct=std_correct;
mean_cellf_err=[nanmean(squeeze(store_f_err(:,1,:)),2),nanmean(squeeze(store_f_err(:,2,:)),2),nanmean(squeeze(store_f_err(:,3,:)),2)];
mean_cellf_err(find(mean_cellf_err==Inf))=NaN;
mean_cellf_err(find(mean_cellf_err>100))=NaN;
main_effect_err=nanmean(mean_cellf_err);
std_err=std(mean_cellf_err,0,'omitnan')/sqrt(size(mean_cellf_err,1));
errhigh_err=std_err;
errlow_err=std_err;
subplot(1,2,1);
bar([1,2,3],main_effect_correct);
hold on;
er = errorbar([1,2,3],main_effect_correct,errlow_correct,errhigh_correct);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlimit = get(gca,'xlim');
title('Partietal location vs epoch,correct');
set(gca,'xticklabel',{'stim location','task epoch','interaction'})
subplot(1,2,2);
bar([1,2,3],main_effect_err);
hold on;
er = errorbar([1,2,3],main_effect_err,errlow_err,errhigh_err);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlimit = get(gca,'xlim');
title('Partietal location vs epoch,err');
set(gca,'xticklabel',{'stim location','task epoch','interaction'})
%{
figure;
subplot(1,3,1);
scatter(mean_cellf_correct(:,1),mean_cellf_err(:,1));
subplot(1,3,2);
scatter(mean_cellf_correct(:,2),mean_cellf_err(:,2));
subplot(1,3,3);
scatter(mean_cellf_correct(:,3),mean_cellf_err(:,3));
%}
figure;
subplot(1,2,1);
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
title('Partietal location vs epoch,correct');
set(gca,'xticklabel',{'stim location','task epoch','interaction'})
%ylim([0,0.3]);
text(1,0.4,strcat('n=',num2str(length(data_inclusion_index))));

subplot(1,2,2);
main_effect_y=[mean(prop_sigstim_err_main1),mean(prop_sigstim_err_main2),mean(prop_sigstim_err_interaction)];
std_y=[std(prop_sigstim_err_main1)/sqrt(50),std(prop_sigstim_err_main2)/sqrt(50),std(prop_sigstim_err_interaction)/sqrt(50)];
errhigh=std_y;
errlow=std_y;
bar([1,2,3],main_effect_y);
hold on;
er = errorbar([1,2,3],main_effect_y,errlow,errhigh);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlimit = get(gca,'xlim');
plot(xlimit,[0.05 0.05]);
title('Parietal location vs epoch,err');
set(gca,'xticklabel',{'stim location','task epoch','interaction'})
%ylim([0,0.3]);
text(1,0.4,strcat('n=',num2str(length(data_inclusion_index))));

