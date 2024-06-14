%find significant cells 
load('all_spatial_info.mat');
load('all_spatial_data.mat');
all_location_names=unique(all_spatial_info(:,5));
group=all_spatial_info(:,4);  % pre post or error trials
%IndexC = strfind(group,'PRE');
IndexC = strfind(group,'POST');
find_index = find(not(cellfun('isempty',IndexC)));
good_data=zeros(1,length(find_index)); %vectore indicate if cell was used

for i=1:length(find_index)

    cell_data=all_spatial_data(find_index(i),:);
    cell_info=all_spatial_info(find_index(i),:);  
    cell_location=cell_info(5);
    cell_location_num=strcmp(cell_location,all_location_names);
    vector_cell_loc(i)=find(cell_location_num==1); %1,Anterior-Dorsal 2, Anterior-Ventral 3, Mid-dorsal 4,Posterior-Dorsal 5,Posterior-Ventral
       
    all_stim=[];
    all_delay=[];
    stim_loc=[];
    task_epoch=[];
    ifmatch=[];
    sample_loc_lookup=[5 6 7 8 1 2 3 4];
    for j=1:8 
        cueclass_data=cell_data{j};
        cue_class_trials(j)=length(cueclass_data);
        if cue_class_trials(j)>0
            if length(fieldnames(cueclass_data))>10
        cueclass_fix_count=[cueclass_data.fix];
        cueclass_cue_rate=[cueclass_data.cuerate];
        cueclass_cuedelay_rate=[cueclass_data.cuedelay];
        cueclass_sample_rate=[cueclass_data.samplerate];
        cueclass_sampledelay_rate=[cueclass_data.sampledelay];
        cueclass_ismatch=[cueclass_data.IsMatch];
            else
        cueclass_fix_count=[];
        cueclass_cue_rate=[];
        cueclass_cuedelay_rate=[];
        cueclass_sample_rate=[];
        cueclass_sampledelay_rate=[];
        cueclass_ismatch=[];
            end
        else
        cueclass_fix_count=[];
        cueclass_cue_rate=[];
        cueclass_cuedelay_rate=[];
        cueclass_sample_rate=[];
        cueclass_sampledelay_rate=[];
        cueclass_ismatch=[];
        end
            
        if length(cueclass_ismatch)==length(cueclass_fix_count)
        %some cells have ifmatch data but no spiking data
        good_data(i)=1;
        cue_loc=j*ones(1,length(cueclass_fix_count));
        ifmatch=[ifmatch,cueclass_ismatch,cueclass_ismatch]; 
        temp_ismatch=cueclass_ismatch;
        sample_loc=cue_loc.*(temp_ismatch)+sample_loc_lookup(cue_loc).*(temp_ismatch*(-1)+1);     
        stim_loc=[stim_loc,cue_loc,sample_loc];
        task_epoch=[task_epoch,[ones(1,length(cueclass_fix_count)),2*ones(1,length(cueclass_fix_count))]];

        all_stim=[all_stim,cueclass_cue_rate,cueclass_sample_rate];
        all_delay=[all_delay,cueclass_cuedelay_rate,cueclass_sampledelay_rate];          
        end
    end
        cell_numtrials(i)=sum(cue_class_trials);
 
    if (length(ifmatch)>0 && isnan(ifmatch(1))~=1)   %some cell have ifmatch column NaN?
    select_matchonly=(ifmatch==1);
    [stim_p(i,:),tbl_stim]=anovan(all_stim(select_matchonly),{stim_loc(select_matchonly) task_epoch(select_matchonly)},'model','interaction','display','off');
    %[stim_p(i,:),tbl_stim]=anovan(all_stim,{stim_loc ifmatch},'model','interaction','display','off');
    f_stim(i,:)=cell2mat(tbl_stim(2:4,6))';
    [delay_p(i,:),tbl_delay]=anovan(all_delay(select_matchonly),{stim_loc(select_matchonly) task_epoch(select_matchonly)},'model','interaction','display','off');
    %[delay_p(i,:),tbl_delay]=anovan(all_delay,{stim_loc ifmatch},'model','interaction','display','off');
    f_delay(i,:)=cell2mat(tbl_delay(2:4,6))';
    else
    good_data(i)=0;
    stim_p(i,:)=NaN;
    delay_p(i,:)=NaN;
    f_stim(i,1:3)=NaN;
    f_delay(i,1:3)=NaN;
    end

    end
    data_inclusion_index=intersect(find(cell_numtrials>95),find(good_data==1));

    find_sigstim_main1=intersect(find(stim_p(:,1)<0.05),data_inclusion_index);
    find_sigstim_main2=intersect(find(stim_p(:,2)<0.05),data_inclusion_index);
    find_sigdelay_main1=intersect(find(delay_p(:,1)<0.05),data_inclusion_index);
    find_sigdelay_main2=intersect(find(delay_p(:,2)<0.05),data_inclusion_index);
    find_sigstim_interaction=intersect(find(stim_p(:,3)<0.05),data_inclusion_index);
    find_sigdelay_interaction=intersect(find(delay_p(:,3)<0.05),data_inclusion_index);


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
       
    find_sigdelay_CS1=setdiff(setdiff(find_sigdelay_main1,find_sigdelay_main2),find_sigdelay_interaction);
    find_sigdelay_CS2=setdiff(setdiff(find_sigdelay_main2,find_sigdelay_main1),find_sigdelay_interaction);
    find_sigdelay_informative=unique([find_sigdelay_main1;find_sigdelay_main2;find_sigdelay_interaction]);
    find_sigdelay_LMS=setdiff(intersect(find_sigdelay_main1,find_sigdelay_main2),find_sigdelay_interaction);
    find_sigdelay_NMS=find_sigdelay_interaction; 
    find_sigdelay_NMS1=intersect(find_sigdelay_interaction,intersect(find_sigdelay_main1,find_sigdelay_main2));
    find_sigdelay_NMS2=intersect(find_sigdelay_interaction,setdiff(find_sigdelay_main1,find_sigdelay_main2));
    find_sigdelay_NMS3=intersect(find_sigdelay_interaction,setdiff(find_sigdelay_main2,find_sigdelay_main1));
    find_sigdelay_NMS4=setdiff(find_sigdelay_interaction,unique([find_sigdelay_main1;find_sigdelay_main2]));
    find_sigdelay_CS=[find_sigdelay_CS1;find_sigdelay_CS2];
    
    prop_sigstim_main1=length(find_sigstim_main1)/length(data_inclusion_index);
    prop_sigstim_main2=length(find_sigstim_main2)/length(data_inclusion_index);
    prop_sigdelay_main1=length(find_sigdelay_main1)/length(data_inclusion_index);
    prop_sigdelay_main2=length(find_sigdelay_main2)/length(data_inclusion_index);
    prop_sigstim_interaction=length(find_sigstim_interaction)/length(data_inclusion_index);
    prop_sigdelay_interaction=length(find_sigdelay_interaction)/length(data_inclusion_index);
    
    find_informative=unique([find_sigstim_main1;find_sigstim_main2;find_sigdelay_main1;find_sigdelay_main2]);
    
    %proportion of informative cells    
    proportion_informative=length(find_informative)/length(data_inclusion_index);
    main_effect_y=[prop_sigstim_main1,prop_sigdelay_main1;prop_sigstim_main2,prop_sigdelay_main2;prop_sigstim_interaction,prop_sigdelay_interaction];
    bar(main_effect_y,'grouped');
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
    title('Pre training');
    %title('Post training');
    set(gca,'xticklabel',{'stim location','task epoch','interaction'})
    legend('stim','delay');
    ylim([0,0.55]);
    text(1,0.4,strcat('n=',num2str(length(data_inclusion_index))));
    
    figure;
    labels={'CS','LMS','NMS'};
    subplot(1,2,1)
    mix_ratio_stim=[length(find_sigstim_CS)/length(find_sigstim_informative),...
        length(find_sigstim_LMS)/length(find_sigstim_informative),...
        length(find_sigstim_NMS)/length(find_sigstim_informative)];
    pie(mix_ratio_stim,labels);
    title(strcat('sample n=',num2str(length(find_sigstim_informative))));
    subplot(1,2,2)
    mix_ratio_delay=[length(find_sigdelay_CS)/length(find_sigdelay_informative),...
        length(find_sigdelay_LMS)/length(find_sigdelay_informative),...
        length(find_sigdelay_NMS)/length(find_sigdelay_informative)];
    pie(mix_ratio_delay,labels);
    title(strcat('delay n=',num2str(length(find_sigdelay_informative))));
    
    figure;
    labels={'CS','LMS','NMS','NS'};
    subplot(1,2,1)
    mix_ratio_stim=[length(find_sigstim_CS)/length(data_inclusion_index),...
        length(find_sigstim_LMS)/length(data_inclusion_index),...
        length(find_sigstim_NMS)/length(data_inclusion_index),...
        (length(data_inclusion_index)-length(find_sigstim_informative))/length(data_inclusion_index)];
    pie(mix_ratio_stim,labels);
    title(strcat('stim n=',num2str(length(data_inclusion_index))));
   
    subplot(1,2,2)
    mix_ratio_delay=[length(find_sigdelay_CS)/length(data_inclusion_index),...
        length(find_sigdelay_LMS)/length(data_inclusion_index),...
        length(find_sigdelay_NMS)/length(data_inclusion_index),...
        (length(data_inclusion_index)-length(find_sigdelay_informative))/length(data_inclusion_index)];
    pie(mix_ratio_delay,labels);
    title(strcat('delay n=',num2str(length(data_inclusion_index))));
%{
     figure;
    subplot(2,2,1);
    hold on;
    scatter(ones(1,length(data_inclusion_index)),f_stim(data_inclusion_index,1), 1);
    scatter(2*ones(1,length(data_inclusion_index)),f_delay(data_inclusion_index,1), 1);
    xticks([1,2]);
    xlim([0,3]);
    xticklabels({'stim','delay'});

  %  ylim([0,100]);
    
    subplot(2,2,3);
    hold on;
    scatter(1,nanmean(f_stim(data_inclusion_index,1)), 80,'x');
    scatter(2,nanmean(f_delay(data_inclusion_index,1)), 80,'x');
    xticks([1,2]);
    xlim([0,3]);
    xticklabels({'stim','delay'});  
    
    subplot(2,2,2);
    hold on;
    scatter(ones(1,length(data_inclusion_index)),f_stim(data_inclusion_index,2), 1);
    scatter(2*ones(1,length(data_inclusion_index)),f_delay(data_inclusion_index,2), 1);
    xticks([1,2]);
    xlim([0,3]);
    xticklabels({'stim','delay'});
  %  ylim([0,60]);
    
    subplot(2,2,4);
    hold on;
    scatter(1,nanmean(f_stim(data_inclusion_index,2)), 80,'x');
    scatter(2,nanmean(f_delay(data_inclusion_index,2)), 80,'x');
    xticks([1,2]);
    xlim([0,3]);
    xticklabels({'stim','delay'});  
%}    
    figure;
plot_area=5;
area_name={'Anterior-Dorsal' 'Anterior-Ventral' 'Mid-dorsal' 'Posterior-Dorsal' 'Posterior-Ventral'};
all_areacells=intersect(find(vector_cell_loc==plot_area),data_inclusion_index);
area_stimmain1=intersect(find_sigstim_main1,find(vector_cell_loc==plot_area));
area_delaymain1=intersect(find_sigdelay_main1,find(vector_cell_loc==plot_area));
area_stimmain2=intersect(find_sigstim_main2,find(vector_cell_loc==plot_area));
area_delaymain2=intersect(find_sigdelay_main2,find(vector_cell_loc==plot_area));
area_stim_interaction=intersect(find_sigstim_interaction,find(vector_cell_loc==plot_area));
area_delay_interaction=intersect(find_sigdelay_interaction,find(vector_cell_loc==plot_area));

prop_area_stimmain1=length(area_stimmain1)/length(all_areacells);
prop_area_delaymain1=length(area_delaymain1)/length(all_areacells);
prop_area_stimmain2=length(area_stimmain2)/length(all_areacells);
prop_area_delaymain2=length(area_delaymain2)/length(all_areacells);
prop_area_stiminter=length(area_stim_interaction)/length(all_areacells);
prop_area_delayinter=length(area_delay_interaction)/length(all_areacells);

main_effect_area=[prop_area_stimmain1,prop_area_delaymain1;prop_area_stimmain2,prop_area_delaymain2;prop_area_stiminter,prop_area_delayinter];
bar(main_effect_area,'grouped');
hold on;
xlimit = get(gca,'xlim');
plot(xlimit,[0.05 0.05]);
ylim([0,0.55]);
title(strcat(area_name{plot_area},' Pre training n=',num2str(length(all_areacells))));
%title(strcat(area_name{plot_area},' Post training n=',num2str(length(all_areacells))));

%set(gca,'xticklabel',{'stim location','epoch','interaction'});
set(gca,'xticklabel',{'stim location','match/nonmatch'});
legend('stim','delay');