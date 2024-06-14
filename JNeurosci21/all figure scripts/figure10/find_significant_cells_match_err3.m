%this version report f score rather than significant proportion
load('all_spatial_info.mat');
load('all_spatial_data.mat');
load('all_spatialerr_info.mat');
load('all_spatialerr_data.mat');
num_class_threshold=[2,3,4,5,6];
num_trial_threshold=[2,3,4,5,6];
post_filter_count=zeros(5,5);
sample_loc_lookup=[5 6 7 8 1 2 3 4];

%check how many neurons pass the requirement
for i=1:length(num_class_threshold)
    %disp(i);
    for j=1:length(num_trial_threshold)
     %   disp(j);
        store_index=[];
        temp_numclass_line=num_class_threshold(i);
        temp_numtrial_line=num_trial_threshold(j);
        for n=1:size(all_spatialerr_data,1)
        cell_correct_data=all_spatial_data(1735+n,:);
        cell_err_data=all_spatialerr_data(n,:);
        for p=1:8
            class_err_data=cell_err_data{p};
            class_correct_data=cell_correct_data{p};
            if (length(class_err_data)>0 & isfield(class_err_data,'IsMatch'))
            err_ismatch=[class_err_data.IsMatch];
            cell_err_nummatch(p)=length(find(err_ismatch));  %in certain class number of nummatch in err trials
            cell_err_numnonmatch(p)=length(err_ismatch)-cell_err_nummatch(p); %in certain class class number of match in err trials
            else
            cell_err_nummatch(p)=0;
            cell_err_numnonmatch(p)=0;
            end
            correct_ismatch=[class_correct_data.IsMatch]; %in certn class number of nummatch in correct trials
            cell_correct_nummatch(p)=length(find(correct_ismatch)); %in certain class number of match in correct trials
            cell_correct_numnonmatch(p)=length(correct_ismatch)-cell_correct_nummatch(p);
        end
        err_numtrial_qualify=(cell_err_nummatch>=temp_numtrial_line).*(circshift(cell_err_numnonmatch,4)>=temp_numtrial_line);
        correct_numtrial_qualify=(cell_correct_nummatch>=temp_numtrial_line).*(circshift(cell_correct_numnonmatch,4)>=temp_numtrial_line);

        if sum(err_numtrial_qualify.*correct_numtrial_qualify)>=temp_numclass_line
        post_filter_count(i,j)=post_filter_count(i,j)+1;
        post_filter_find{i,j,n}=find(err_numtrial_qualify.*correct_numtrial_qualify); %store class that qualify for each threshold conditiona and cell
        store_index=[store_index,n];  %store cell index that qualify
        end
        end
        post_filter_index{i,j}=store_index;
    end
end

select_class_threshold=2;
select_trial_threshold=2; 
after_filter_index=post_filter_index{select_class_threshold,select_trial_threshold};

find_index = after_filter_index;
class_to_use=squeeze(post_filter_find(select_class_threshold,select_trial_threshold,:));
for n=1:50
disp(n);
for i=1:length(find_index)

     %number of resample
    cell_class_touse=class_to_use{find_index(i)};
    cell_data=all_spatial_data(1735+find_index(i),:);
    cell_info=all_spatial_info(1735+find_index(i),:); 
    cell_err_data=all_spatialerr_data(find_index(i),:);


       
    %all_fix=[];
    %all_cue=[];
    %all_cuedelay=[];
    all_sample=[];
    all_sampledelay=[];
    
    %all_errfix=[];
    %all_errcue=[];
    %all_errcuedelay=[];
    all_errsample=[];
    all_errsampledelay=[];
    
 %   all_cueloc=[];
    all_sampleloc=[];
    all_ifmatch=[];

    for j=1:length(cell_class_touse)
        
        cueclass_err_data1=cell_err_data{cell_class_touse(j)};        
        cueclass_err_data2=cell_err_data{sample_loc_lookup(cell_class_touse(j))};
        cueclass_data1=cell_data{cell_class_touse(j)};
        cueclass_data2=cell_data{sample_loc_lookup(cell_class_touse(j))};

        
            err_ismatch1=[cueclass_err_data1.IsMatch];% ismatch data for err trials
            err_ismatch2=[cueclass_err_data2.IsMatch];
            correct_ismatch1=[cueclass_data1.IsMatch];  %ismatch data for correct trials
            correct_ismatch2=[cueclass_data2.IsMatch];
            
            class_err_nummatch=length(find(err_ismatch1));
            class_err_findmatch=find(err_ismatch1);
            class_err_numnonmatch=length(find(~err_ismatch2));
            class_err_findnonmatch=find(~err_ismatch2);            
            
            class_nummatch=length(find(correct_ismatch1));
            class_findmatch=find(correct_ismatch1);
            class_numnonmatch=length(find(~correct_ismatch2));
            class_findnonmatch=find(~correct_ismatch2);
            
            lower_bound=min([class_nummatch,class_err_nummatch,class_numnonmatch,class_err_numnonmatch]);
            
             temp_rand_select=randperm(class_err_nummatch); 
             temp_errmatch_select=temp_rand_select(1:lower_bound);
             temp_rand_select=randperm(class_err_numnonmatch); 
             temp_errnonmatch_select=temp_rand_select(1:lower_bound);
             temp_rand_select=randperm(class_nummatch); 
             temp_match_select=temp_rand_select(1:lower_bound);
             temp_rand_select=randperm(class_numnonmatch); 
             temp_nonmatch_select=temp_rand_select(1:lower_bound);
             
         %    cueclass_errfix_count=[cueclass_err_data.fix]; %fix count for err trials
         %    cueclass_errfix_count=[cueclass_errfix_count(class_err_findmatch(temp_errmatch_select)),cueclass_errfix_count(class_err_findnonmatch(temp_errnonmatch_select))];  
         %    cueclass_errcue_count=[cueclass_err_data.cuerate];
         %    cueclass_errcue_count=[cueclass_errcue_count(class_err_findmatch(temp_errmatch_select)),cueclass_errcue_count(class_err_findnonmatch(temp_errnonmatch_select))];
         %    cueclass_errcuedelay_count=[cueclass_err_data.cuedelay];
         %    cueclass_errcuedelay_count=[cueclass_errcuedelay_count(class_err_findmatch(temp_errmatch_select)),cueclass_errcuedelay_count(class_err_findnonmatch(temp_errnonmatch_select))];
             cueclass_errsample_count1=[cueclass_err_data1.samplerate];
             cueclass_errsample_count2=[cueclass_err_data2.samplerate];
             cueclass_errsample_count=[cueclass_errsample_count1(class_err_findmatch(temp_errmatch_select)),cueclass_errsample_count2(class_err_findnonmatch(temp_errnonmatch_select))];
             cueclass_errsampledelay_count1=[cueclass_err_data1.sampledelay];
             cueclass_errsampledelay_count2=[cueclass_err_data2.sampledelay];
             cueclass_errsampledelay_count=[cueclass_errsampledelay_count1(class_err_findmatch(temp_errmatch_select)),cueclass_errsampledelay_count2(class_err_findnonmatch(temp_errnonmatch_select))];
             cueclass_errismatch1=[cueclass_err_data1.IsMatch];
             cueclass_errismatch2=[cueclass_err_data2.IsMatch];
             cueclass_errismatch=[cueclass_errismatch1(class_err_findmatch(temp_errmatch_select)),cueclass_errismatch2(class_err_findnonmatch(temp_errnonmatch_select))];

         
             %picking trials from correct data

          %   cueclass_fix_count=[cueclass_data1.fix];
          %   cueclass_fix_count=[cueclass_fix_count(class_findmatch(temp_match_select)),cueclass_fix_count(class_findnonmatch(temp_nonmatch_select))];  
          %   cueclass_cue_count=[cueclass_data.cuerate];
          %   cueclass_cue_count=[cueclass_cue_count(class_findmatch(temp_match_select)),cueclass_cue_count(class_findnonmatch(temp_nonmatch_select))];  
          %   cueclass_cuedelay_count=[cueclass_data.cuedelay];
          %   cueclass_cuedelay_count=[cueclass_cuedelay_count(class_findmatch(temp_match_select)),cueclass_cuedelay_count(class_findnonmatch(temp_nonmatch_select))];  
             cueclass_sample_count1=[cueclass_data1.samplerate];
             cueclass_sample_count2=[cueclass_data2.samplerate];
             cueclass_sample_count=[cueclass_sample_count1(class_findmatch(temp_match_select)),cueclass_sample_count2(class_findnonmatch(temp_nonmatch_select))];  
             cueclass_sampledelay_count1=[cueclass_data1.sampledelay];
             cueclass_sampledelay_count2=[cueclass_data2.sampledelay];
             cueclass_sampledelay_count=[cueclass_sampledelay_count1(class_findmatch(temp_match_select)),cueclass_sampledelay_count2(class_findnonmatch(temp_nonmatch_select))];  
             cueclass_ismatch1=[cueclass_data1.IsMatch];
             cueclass_ismatch2=[cueclass_data2.IsMatch];
             cueclass_ismatch=[cueclass_ismatch1(class_findmatch(temp_match_select)),cueclass_ismatch2(class_findnonmatch(temp_nonmatch_select))];
               
            
    %    cue_loc=j*ones(1,lower_bound*2);
    %    all_cueloc=[all_cueloc,cue_loc];
        all_ifmatch=[all_ifmatch,cueclass_errismatch]; 
    %    temp_ismatch=cueclass_errismatch;
    %    sample_loc=cue_loc.*(temp_ismatch)+sample_loc_lookup(cue_loc).*(temp_ismatch*(-1)+1);     
        sample_loc=[j*ones(1,lower_bound*2)];
        all_sampleloc=[all_sampleloc,sample_loc];
        
    %    all_fix=[all_fix,cueclass_fix_count];
    %    all_cue=[all_cue,cueclass_cue_count];
    %    all_cuedelay=[all_cuedelay,cueclass_cuedelay_count];
        all_sample=[all_sample,cueclass_sample_count];
        all_sampledelay=[all_sampledelay,cueclass_sampledelay_count];
  
    %    all_errfix=[all_errfix,cueclass_errfix_count];
    %    all_errcue=[all_errcue,cueclass_errcue_count];
    %    all_errcuedelay=[all_errcuedelay,cueclass_errcuedelay_count];
        all_errsample=[all_errsample,cueclass_errsample_count];
        all_errsampledelay=[all_errsampledelay,cueclass_errsampledelay_count]; 
    end
 
 
   % if length(all_stim)>0 
  %  [fix_p(i,:),tbl_fix]=anovan(all_fix,{all_cueloc all_ifmatch},'model','interaction','display','off');
  %  [errfix_p(i,:),tbl_errfix]=anovan(all_errfix,{all_cueloc all_ifmatch},'model','interaction','display','off');
  %  [cue_p(i,:),tbl_cue]=anovan(all_cue,{all_cueloc all_ifmatch},'model','interaction','display','off');
  %  [errcue_p(i,:),tbl_errcue]=anovan(all_errcue,{all_cueloc all_ifmatch},'model','interaction','display','off');
  %  [cuedelay_p(i,:),tbl_cuedelay]=anovan(all_cuedelay,{all_cueloc all_ifmatch},'model','interaction','display','off');
  %  [errcuedelay_p(i,:),tbl_errcuedelay]=anovan(all_errcuedelay,{all_cueloc all_ifmatch},'model','interaction','display','off');
    [sample_p(i,:),tbl_sample]=anovan(all_sample,{all_sampleloc all_ifmatch},'model','interaction','display','off');
    f_sample(i,:)=cell2mat(tbl_sample(2:4,6))';
    [errsample_p(i,:),tbl_errsample]=anovan(all_errsample,{all_sampleloc all_ifmatch},'model','interaction','display','off');
    f_errsample(i,:)=cell2mat(tbl_errsample(2:4,6))';
    [sampledelay_p(i,:),tbl_sampledelay]=anovan(all_sampledelay,{all_sampleloc all_ifmatch},'model','interaction','display','off');
    f_sampledelay(i,:)=cell2mat(tbl_sampledelay(2:4,6))';
    [errsampledelay_p(i,:),tbl_errsampledelay]=anovan(all_errsampledelay,{all_sampleloc all_ifmatch},'model','interaction','display','off');
    f_errsampledelay(i,:)=cell2mat(tbl_errsampledelay(2:4,6))';
end
  %  store_resample{n,1}=fix_p;
  %  store_resample{n,2}=cue_p;
  %  store_resample{n,3}=cuedelay_p;
    store_resample{n,1}=f_sample;
    store_resample{n,2}=f_sampledelay;
  %  store_resample{n,6}=errfix_p;
  %  store_resample{n,7}=errcue_p;
  %  store_resample{n,8}=errcuedelay_p;
    store_resample{n,3}=f_errsample;
    store_resample{n,4}=f_errsampledelay;
end
for n=1:50
    sample_f=store_resample{n,1};
    sampledelay_f=store_resample{n,2};
    errsample_f=store_resample{n,3};
    errsampledelay_f=store_resample{n,4};

    mean_sample_main1(n,:)=sample_f(:,1);
    mean_sample_main2(n,:)=sample_f(:,2);
    mean_sample_interaction(n,:)=sample_f(:,3);
    mean_sampledelay_main1(n,:)=sampledelay_f(:,1);
    mean_sampledelay_main2(n,:)=sampledelay_f(:,2);
    mean_sampledelay_interaction(n,:)=sampledelay_f(:,3);
    
    mean_errsample_main1(n,:)=errsample_f(:,1);
    mean_errsample_main2(n,:)=errsample_f(:,2);
    mean_errsample_interaction(n,:)=errsample_f(:,3);
    mean_errsampledelay_main1(n,:)=errsampledelay_f(:,1);
    mean_errsampledelay_main2(n,:)=errsampledelay_f(:,2);
    mean_errsampledelay_interaction(n,:)=errsampledelay_f(:,3);

end
figure;
  %  main_effect_y=[mean(prop_sigfix_main1),mean(prop_sigcue_main1),mean(prop_sigcuedelay_main1),mean(prop_sigsample_main1),mean(prop_sigsampledelay_main1);...
  %      mean(prop_sigfix_main2),mean(prop_sigcue_main2),mean(prop_sigcuedelay_main2),mean(prop_sigsample_main2),mean(prop_sigsampledelay_main2);...
  %      mean(prop_sigfix_interaction),mean(prop_sigcue_interaction),mean(prop_sigcuedelay_interaction),mean(prop_sigsample_interaction),mean(prop_sigsampledelay_interaction)];
     main_effect_y=[mean(mean_sample_main1),mean(mean_sampledelay_main1);...
        mean(mean_sample_main2),mean(mean_sampledelay_main2);...
        mean(mean_sample_interaction),mean(mean_sampledelay_interaction)];  
  subplot(1,2,1);
    bar(main_effect_y,'grouped');
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
    title('correct dataset');
    set(gca,'xticklabel',{'stim location','match/nonmatch','interaction'})
    legend('stim','delay');
    ylim([0,2]);
 %   text(1,0.2,strcat('n=',num2str(length(data_inclusion_index))));
    
 %   main_effect_y=[mean(prop_sigerrfix_main1),mean(prop_sigerrcue_main1),mean(prop_sigerrcuedelay_main1),mean(prop_sigerrsample_main1),mean(prop_sigerrsampledelay_main1);...
 %       mean(prop_sigerrfix_main2),mean(prop_sigerrcue_main2),mean(prop_sigerrcuedelay_main2),mean(prop_sigerrsample_main2),mean(prop_sigerrsampledelay_main2);...
 %       mean(prop_sigerrfix_interaction),mean(prop_sigerrcue_interaction),mean(prop_sigerrcuedelay_interaction),mean(prop_sigerrsample_interaction),mean(prop_sigerrsampledelay_interaction)];
     main_effect_y=[mean(mean_errsample_main1),mean(mean_errsampledelay_main1);...
        mean(mean_errsample_main2),mean(mean_errsampledelay_main2);...
        mean(mean_errsample_interaction),mean(mean_errsampledelay_interaction)];   
    subplot(1,2,2);
    bar(main_effect_y,'grouped');
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
    title('error dataset');
    set(gca,'xticklabel',{'stim location','match/nonmatch','interaction'})
    legend('stim','delay');
    ylim([0,2]);
 %   text(1,0.2,strcat('n=',num2str(length(data_inclusion_index))));
    