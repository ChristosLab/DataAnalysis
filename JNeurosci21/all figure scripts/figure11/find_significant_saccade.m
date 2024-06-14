%find significant cells 

%load('overlap_spatial_info.mat');
%load('overlap_spatial_data.mat');
%load('3task_spatial.mat');  %analyse only 3 task cells
%overlap_spatial_data=overlap_spatial_data(spatial_96touse,:);
%overlap_spatial_info=overlap_spatial_info(spatial_96touse,:);

%all_location_names=unique(overlap_spatial_info(:,5));
%group=overlap_spatial_info(:,4);  
load('all_spatial_info.mat');
load('all_spatial_data.mat');
load('spatial_sacdirect.mat');
all_location_names=unique(all_spatial_info(:,5));
%extract animal ID
all_cellname=all_spatial_info(:,1);
for i=1:length(all_cellname)
    temp_name=all_cellname{i};
    all_animalID{i}=temp_name(1:3);
end
unique_animalID=unique(lower(all_animalID));
IndexA=strfind(lower(all_animalID),'nin'); 
group=all_spatial_info(:,4); 
%IndexC = strfind(group,'PRE');
IndexC = strfind(group,'POST');%only post trainig have saccade
%find_index = intersect(find(not(cellfun('isempty',IndexC))),find(not(cellfun('isempty',IndexA)))); %only use cells from certain animal 
find_index = find(not(cellfun('isempty',IndexC))); %use all pre or post cells
good_data=zeros(1,length(find_index)); %vectore indicate if cell was used
edge1=[-1,0,0.5,2,2.5,4.2];   
edge2=[0,0.5,2,2.5,4,5];
%skip_index=[2860:2865]; %behavior and spike data does not match
for i=1:length(find_index)
     cell_data=all_spatial_data(find_index(i),:);
     cell_info=all_spatial_info(find_index(i),:);
     cell_sacdirect=all_trialdirect{find_index(i)};
   % cell_data=overlap_spatial_data(find_index(i),:);
   % cell_info=overlap_spatial_info(find_index(i),:);
    cell_location=cell_info(5);
    cell_location_num=strcmp(cell_location,all_location_names);
    vector_cell_loc(i)=find(cell_location_num==1); %1,Anterior-Dorsal 2, Anterior-Ventral 3, Mid-dorsal 4,Posterior-Dorsal 5,Posterior-Ventral
    
    all_fix=[];
    all_cue=[];
    all_cuedelay=[];
    all_sample=[];
    all_sampledelay=[];
    all_choice=[];
    sac_loc=[];
    ifmatch=[];
    cue_class_trials=[];
    if i==27
       disp('test');
    end
    for j=1:8        
        cueclass_data=cell_data{j};
        cue_class_trials(j)=length(cueclass_data);
        if cue_class_trials(j)>0
            if length(fieldnames(cueclass_data))>10
        cueclass_trialnum=[cueclass_data.trialnum];
        if max(cueclass_trialnum)<=length(cell_sacdirect)
        cueclass_sacdirect=cell_sacdirect(cueclass_trialnum);
        else
            break;
        end
        cueclass_ismatch=[cueclass_data.IsMatch];
        cueclass_spiketimes={cueclass_data.TS};
        cueclass_ontime=[cueclass_data.Cue_onT];
        cueclass_bincount=[];
               for p=1:cue_class_trials(j) %loop through correct trials
                 temp_TS=cueclass_spiketimes{p}-cueclass_ontime(p);
                 for q=1:length(edge1)  %loop through timebin            
                   cueclass_bincount(p,q)=length(find(temp_TS>edge1(q) & temp_TS<edge2(q)));
                 end
               end
        cueclass_fix_count=cueclass_bincount(:,1)';
        cueclass_cue_rate=cueclass_bincount(:,2)';
        cueclass_cuedelay_rate=cueclass_bincount(:,3)';
        cueclass_sample_rate=cueclass_bincount(:,4)';
        cueclass_sampledelay_rate=cueclass_bincount(:,5)';
        cueclass_choice_rate=cueclass_bincount(:,6)';
            else
        cueclass_fix_count=[];
        cueclass_cue_rate=[];
        cueclass_cuedelay_rate=[];
        cueclass_sample_rate=[];
        cueclass_sampledelay_rate=[];
        cueclass_choice_rate=[];
        cueclass_ismatch=[];
            end
        else
        cueclass_fix_count=[];
        cueclass_cue_rate=[];
        cueclass_cuedelay_rate=[];
        cueclass_sample_rate=[];
        cueclass_sampledelay_rate=[];
        cueclass_choice_rate=[];
        cueclass_ismatch=[];
        end
        %two-way anova 
        if length(cueclass_ismatch)==length(cueclass_fix_count)
        %found some cells have ifmatch data but no spiking data
        good_data(i)=1;
        ifmatch=[ifmatch,cueclass_ismatch]; 
        sac_loc=[sac_loc,cueclass_sacdirect]; 
        all_fix=[all_fix,cueclass_fix_count];
        all_cue=[all_cue,cueclass_cue_rate];
        all_cuedelay=[all_cuedelay,cueclass_cuedelay_rate];
        all_sample=[all_sample,cueclass_sample_rate];
        all_sampledelay=[all_sampledelay,cueclass_sampledelay_rate]; 
        all_choice=[all_choice,cueclass_choice_rate];
        end
    end
        cell_numtrials(i)=sum(cue_class_trials);

    % two-way anova for interaction, it is redundant with both one way and
    % twoway just to cross check results
    if (length(ifmatch)>0 && isnan(ifmatch(1))~=1)  %why some cell have ifmatch column NaN?
    [interaction_fix_p(i,:),tbl_fix]=anovan(all_fix,{sac_loc ifmatch},'model','interaction','display','off');
    f_fix(i,:)=cell2mat(tbl_fix(2:4,6))';
    [interaction_cue_p(i,:),tbl_cue]=anovan(all_cue,{sac_loc ifmatch},'model','interaction','display','off');
    f_cue(i,:)=cell2mat(tbl_cue(2:4,6))';
    [interaction_cuedelay_p(i,:),tbl_cuedelay]=anovan(all_cuedelay,{sac_loc ifmatch},'model','interaction','display','off');
    f_cuedelay(i,:)=cell2mat(tbl_cuedelay(2:4,6))';
    [interaction_sample_p(i,:),tbl_sample]=anovan(all_sample,{sac_loc ifmatch},'model','interaction','display','off');
    f_sample(i,:)=cell2mat(tbl_sample(2:4,6))';
    [interaction_sampledelay_p(i,:),tbl_sampledelay]=anovan(all_sampledelay,{sac_loc ifmatch},'model','interaction','display','off');
    f_sampledelay(i,:)=cell2mat(tbl_sampledelay(2:4,6))';
    [interaction_choice_p(i,:),tbl_choice]=anovan(all_choice,{sac_loc ifmatch},'model','interaction','display','off');
    f_choice(i,:)=cell2mat(tbl_choice(2:4,6))';
    else
    good_data(i)=0;
    interaction_fix_p(i,1:3)=NaN;
    interaction_cue_p(i,1:3)=NaN;
    interaction_cuedelay_p(i,1:3)=NaN;
    interaction_sample_p(i,1:3)=NaN;
    interaction_sampledelay_p(i,1:3)=NaN;
    f_fix(i,1:3)=NaN;
    f_cue(i,1:3)=NaN;
    f_cuedelay(i,1:3)=NaN;
    f_sample(i,1:3)=NaN;
    f_sampledelay(i,1:3)=NaN;
    end

end
    data_inclusion_index=intersect(find(cell_numtrials>95),find(good_data==1));

    find_sigfix_main1=intersect(find(interaction_fix_p(:,1)<0.05),data_inclusion_index);
    find_sigcue_main1=intersect(find(interaction_cue_p(:,1)<0.05),data_inclusion_index);
    find_sigcuedelay_main1=intersect(find(interaction_cuedelay_p(:,1)<0.05),data_inclusion_index);
    find_sigsample_main1=intersect(find(interaction_sample_p(:,1)<0.05),data_inclusion_index);
    find_sigsampledelay_main1=intersect(find(interaction_sampledelay_p(:,1)<0.05),data_inclusion_index);
    find_sigchoice_main1=intersect(find(interaction_choice_p(:,1)<0.05),data_inclusion_index);
      
    find_sigfix_main2=intersect(find(interaction_fix_p(:,2)<0.05),data_inclusion_index);
    find_sigcue_main2=intersect(find(interaction_cue_p(:,2)<0.05),data_inclusion_index);
    find_sigcuedelay_main2=intersect(find(interaction_cuedelay_p(:,2)<0.05),data_inclusion_index);
    find_sigsample_main2=intersect(find(interaction_sample_p(:,2)<0.05),data_inclusion_index);
    find_sigsampledelay_main2=intersect(find(interaction_sampledelay_p(:,2)<0.05),data_inclusion_index);
    find_sigchoice_main2=intersect(find(interaction_choice_p(:,2)<0.05),data_inclusion_index);
        
    find_sigfix_interaction=intersect(find(interaction_fix_p(:,3)<0.05),data_inclusion_index);
    find_sigcue_interaction=intersect(find(interaction_cue_p(:,3)<0.05),data_inclusion_index);
    find_sigcuedelay_interaction=intersect(find(interaction_cuedelay_p(:,3)<0.05),data_inclusion_index);
    find_sigsample_interaction=intersect(find(interaction_sample_p(:,3)<0.05),data_inclusion_index);
    find_sigsampledelay_interaction=intersect(find(interaction_sampledelay_p(:,3)<0.05),data_inclusion_index);
    find_sigchoice_interaction=intersect(find(interaction_choice_p(:,3)<0.05),data_inclusion_index);

    
    find_sigchoice_CS1=setdiff(setdiff(find_sigchoice_main1,find_sigchoice_main2),find_sigchoice_interaction);
    find_sigchoice_CS2=setdiff(setdiff(find_sigchoice_main2,find_sigchoice_main1),find_sigchoice_interaction);
    find_sigchoice_informative=unique([find_sigchoice_main1;find_sigchoice_main2;find_sigchoice_interaction]);
    find_sigchoice_LMS=setdiff(intersect(find_sigchoice_main1,find_sigchoice_main2),find_sigchoice_interaction);
    find_sigchoice_NMS=find_sigchoice_interaction;
    find_sigchoice_NMS1=intersect(find_sigchoice_interaction,intersect(find_sigchoice_main1,find_sigchoice_main2));
    find_sigchoice_NMS2=intersect(find_sigchoice_interaction,setdiff(find_sigchoice_main1,find_sigchoice_main2));
    find_sigchoice_NMS3=intersect(find_sigchoice_interaction,setdiff(find_sigchoice_main2,find_sigchoice_main1));
    find_sigchoice_NMS4=setdiff(find_sigchoice_interaction,unique([find_sigchoice_main1;find_sigchoice_main2]));
    find_sigchoice_CS=[find_sigchoice_CS1;find_sigchoice_CS2];
    
    
    prop_sigfix_main1=length(find_sigfix_main1)/length(data_inclusion_index);
    prop_sigcue_main1=length(find_sigcue_main1)/length(data_inclusion_index);
    prop_sigcuedelay_main1=length(find_sigcuedelay_main1)/length(data_inclusion_index);
    prop_sigsample_main1=length(find_sigsample_main1)/length(data_inclusion_index);
    prop_sigsampledelay_main1=length(find_sigsampledelay_main1)/length(data_inclusion_index);
    prop_sigchoice_main1=length(find_sigchoice_main1)/length(data_inclusion_index);
    
    prop_sigfix_main2=length(find_sigfix_main2)/length(data_inclusion_index);
    prop_sigcue_main2=length(find_sigcue_main2)/length(data_inclusion_index);
    prop_sigcuedelay_main2=length(find_sigcuedelay_main2)/length(data_inclusion_index);
    prop_sigsample_main2=length(find_sigsample_main2)/length(data_inclusion_index);
    prop_sigsampledelay_main2=length(find_sigsampledelay_main2)/length(data_inclusion_index);
    prop_sigchoice_main2=length(find_sigchoice_main2)/length(data_inclusion_index);

    
    prop_sigfix_interaction=length(find_sigfix_interaction)/length(data_inclusion_index);
    prop_sigcue_interaction=length(find_sigcue_interaction)/length(data_inclusion_index);
    prop_sigcuedelay_interaction=length(find_sigcuedelay_interaction)/length(data_inclusion_index);
    prop_sigsample_interaction=length(find_sigsample_interaction)/length(data_inclusion_index);
    prop_sigsampledelay_interaction=length(find_sigsampledelay_interaction)/length(data_inclusion_index);
    prop_sigchoice_interaction=length(find_sigchoice_interaction)/length(data_inclusion_index);

    
    figure;
    main_effect_y=[prop_sigfix_main1,prop_sigcue_main1,prop_sigcuedelay_main1,prop_sigsample_main1,prop_sigsampledelay_main1,prop_sigchoice_main1;...
        prop_sigfix_main2,prop_sigcue_main2,prop_sigcuedelay_main2,prop_sigsample_main2,prop_sigsampledelay_main2,prop_sigchoice_main2;...
        prop_sigfix_interaction,prop_sigcue_interaction,prop_sigcuedelay_interaction,prop_sigsample_interaction,prop_sigsampledelay_interaction,prop_sigchoice_interaction];
    bar(main_effect_y,'grouped');
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
    ylim([0,0.25]);
    text(1,0.2,strcat('n=',num2str(length(data_inclusion_index))));
    %title('Pre training');
    title('Post training');
    set(gca,'xticklabel',{'saccade target location','match/non-match','interaction'})
    legend('fix','cue','cuedelay','sample','sampledelay','choice');
 