%this version report f score rather than significant proportion

load('all_spatial_info.mat');
load('all_spatial_data.mat');
load('all_spatialerr_info.mat');
load('all_spatialerr_data.mat');
num_class_threshold=[2,3,4,5,6];
num_trial_threshold=[2,3,4,5,6];
post_filter_count=zeros(5,5);
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
            cell_err_nummatch(p)=length(find(err_ismatch));
            else
            cell_err_nummatch(p)=0;
            end
            correct_ismatch=[class_correct_data.IsMatch];
            cell_correct_nummatch(p)=length(find(correct_ismatch));
        end

        if sum((cell_err_nummatch>=temp_numtrial_line).*(cell_correct_nummatch>=temp_numtrial_line))>=temp_numclass_line
        post_filter_count(i,j)=post_filter_count(i,j)+1;
        post_filter_find{i,j,n}=find((cell_err_nummatch>=temp_numtrial_line).*(cell_correct_nummatch>=temp_numtrial_line));
        store_index=[store_index,n];
        end
        end
        post_filter_index{i,j}=store_index;
    end
end
%post filter index stores the cells that pass the filter
%post filter find stores the class number to use for past filter cells
select_class_threshold=3;
select_trial_threshold=3; 
after_filter_index=post_filter_index{select_class_threshold,select_trial_threshold};

find_index = after_filter_index;
class_to_use=squeeze(post_filter_find(select_class_threshold,select_trial_threshold,:));
for n=1:50  %resample process
    disp(n);
for i=1:length(find_index)
     %number of resample
    cell_class_touse=class_to_use{find_index(i)};
    cell_data=all_spatial_data(1735+find_index(i),:);
    cell_info=all_spatial_info(1735+find_index(i),:); 
    cell_err_data=all_spatialerr_data(find_index(i),:);

  %  cell_location=cell_info(5);
  %  cell_location_num=strcmp(cell_location,all_location_names);
  %  vector_cell_loc(i)=find(cell_location_num==1); %1,Anterior-Dorsal 2, Anterior-Ventral 3, Mid-dorsal 4,Posterior-Dorsal 5,Posterior-Ventral
       
    all_errstim=[];
    all_errdelay=[];
    all_stim=[];
    all_delay=[];
    stim_loc=[];
    task_epoch=[];
    ifmatch=[];
    sample_loc_lookup=[5 6 7 8 1 2 3 4];

    for j=1:length(cell_class_touse)
        %load correct and error dataset
        cueclass_err_data=cell_err_data{cell_class_touse(j)};
        cueclass_data=cell_data{cell_class_touse(j)};

        %calculate number of trials to pick for certain class
            err_ismatch=[cueclass_err_data.IsMatch];
            correct_ismatch=[cueclass_data.IsMatch];
            class_err_nummatch=length(find(err_ismatch));
            class_err_findmatch=find(err_ismatch);
            class_nummatch=length(find(correct_ismatch));
            class_findmatch=find(correct_ismatch);
            lower_bound=min([class_nummatch,class_err_nummatch]);
            
        %  random select trial number of lower_bound from both dataset
             temp_rand_select=randperm(class_err_nummatch); 
             temp_rand_select=temp_rand_select(1:lower_bound);
             cueclass_errfix_count=[cueclass_err_data.fix];
             cueclass_errfix_count=cueclass_errfix_count(class_err_findmatch(temp_rand_select));
             cueclass_errcue_rate=[cueclass_err_data.cuerate];
             cueclass_errcue_rate=cueclass_errcue_rate(class_err_findmatch(temp_rand_select));
             cueclass_errcuedelay_rate=[cueclass_err_data.cuedelay];
             cueclass_errcuedelay_rate=cueclass_errcuedelay_rate(class_err_findmatch(temp_rand_select));
             cueclass_errsample_rate=[cueclass_err_data.samplerate];
             cueclass_errsample_rate=cueclass_errsample_rate(class_err_findmatch(temp_rand_select));
             cueclass_errsampledelay_rate=[cueclass_err_data.sampledelay];
             cueclass_errsampledelay_rate=cueclass_errsampledelay_rate(class_err_findmatch(temp_rand_select));
             cueclass_errismatch=[cueclass_err_data.IsMatch];
             cueclass_errismatch=cueclass_errismatch(class_err_findmatch(temp_rand_select));
 
             %picking trials from correct data
             temp_rand_select=randperm(class_nummatch);
             temp_rand_select=temp_rand_select(1:lower_bound);
             cueclass_fix_count=[cueclass_data.fix];
             cueclass_fix_count=cueclass_fix_count(class_findmatch(temp_rand_select));
             cueclass_cue_rate=[cueclass_data.cuerate];
             cueclass_cue_rate=cueclass_cue_rate(class_findmatch(temp_rand_select));
             cueclass_cuedelay_rate=[cueclass_data.cuedelay];
             cueclass_cuedelay_rate=cueclass_cuedelay_rate(class_findmatch(temp_rand_select));
             cueclass_sample_rate=[cueclass_data.samplerate];
             cueclass_sample_rate=cueclass_sample_rate(class_findmatch(temp_rand_select));
             cueclass_sampledelay_rate=[cueclass_data.sampledelay];
             cueclass_sampledelay_rate=cueclass_sampledelay_rate(class_findmatch(temp_rand_select));
             cueclass_ismatch=[cueclass_data.IsMatch];
             cueclass_ismatch=cueclass_ismatch(class_findmatch(temp_rand_select));
               
            
        cue_loc=j*ones(1,length(cueclass_errfix_count));
        ifmatch=[ifmatch,cueclass_errismatch,cueclass_errismatch]; 
        temp_ismatch=cueclass_errismatch;
        sample_loc=cue_loc.*(temp_ismatch)+sample_loc_lookup(cue_loc).*(temp_ismatch*(-1)+1);     
        stim_loc=[stim_loc,cue_loc,sample_loc];
        task_epoch=[task_epoch,[ones(1,length(cueclass_errfix_count)),2*ones(1,length(cueclass_errfix_count))]];

        all_errstim=[all_errstim,cueclass_errcue_rate,cueclass_errsample_rate];
        all_errdelay=[all_errdelay,cueclass_errcuedelay_rate,cueclass_errsampledelay_rate]; 
        all_stim=[all_stim,cueclass_cue_rate,cueclass_sample_rate];
        all_delay=[all_delay,cueclass_cuedelay_rate,cueclass_sampledelay_rate];    
  
    end

   % if length(all_stim)>0 
    [stim_p(i,:),tbl_stim]=anovan(all_stim,{stim_loc task_epoch},'model','interaction','display','off');
    f_stim(i,:)=cell2mat(tbl_stim(2:4,6))';
    [errstim_p(i,:),tbl_errstim]=anovan(all_errstim,{stim_loc task_epoch},'model','interaction','display','off');
    f_errstim(i,:)=cell2mat(tbl_errstim(2:4,6))';
    [delay_p(i,:),tbl_delay]=anovan(all_delay,{stim_loc task_epoch},'model','interaction','display','off');
    f_delay(i,:)=cell2mat(tbl_delay(2:4,6))';
    [errdelay_p(i,:),tbl_errdelay]=anovan(all_errdelay,{stim_loc task_epoch},'model','interaction','display','off');
    f_errdelay(i,:)=cell2mat(tbl_errdelay(2:4,6))';

end
    store_resample{n,1}=f_stim;
    store_resample{n,2}=f_delay;
    store_resample{n,3}=f_errstim;
    store_resample{n,4}=f_errdelay;
end
for n=1:50
    stim_f=store_resample{n,1};
    delay_f=store_resample{n,2};
    errstim_f=store_resample{n,3};
    errdelay_f=store_resample{n,4};

    mean_stim_main1(n,:)=stim_f(:,1);
    mean_stim_main2(n,:)=stim_f(:,2);
    mean_stim_interaction(n,:)=stim_f(:,3);
    mean_delay_main1(n,:)=delay_f(:,1);
    mean_delay_main2(n,:)=delay_f(:,2);
    mean_delay_interaction(n,:)=delay_f(:,3);
    
    mean_errstim_main1(n,:)=errstim_f(:,1);
    mean_errstim_main2(n,:)=errstim_f(:,2);
    mean_errstim_interaction(n,:)=errstim_f(:,3);
    mean_errdelay_main1(n,:)=errdelay_f(:,1);
    mean_errdelay_main2(n,:)=errdelay_f(:,2);
    mean_errdelay_interaction(n,:)=errdelay_f(:,3);

end
figure;
    main_effect_y=[mean(mean_stim_main1),mean(mean_delay_main1);mean(mean_stim_main2),mean(mean_delay_main2);mean(mean_stim_interaction),mean(mean_delay_interaction)];
    subplot(1,2,1);
    bar(main_effect_y,'grouped');
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
    title('correct dataset');
    set(gca,'xticklabel',{'stim location','task epoch','interaction'})
    legend('stim','delay');
    ylim([0,3]);
    %text(1,0.4,strcat('n=',num2str(length(data_inclusion_index))));
    
    main_effect_y=[mean(mean_errstim_main1),mean(mean_errdelay_main1);mean(mean_errstim_main2),mean(mean_errdelay_main2);mean(mean_errstim_interaction),mean(mean_errdelay_interaction)];
    subplot(1,2,2);
    bar(main_effect_y,'grouped');
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
    title('error dataset');
    set(gca,'xticklabel',{'stim location','task epoch','interaction'})
    legend('stim','delay');
    ylim([0,3]);
  %  text(1,0.4,strcat('n=',num2str(length(data_inclusion_index))));
    