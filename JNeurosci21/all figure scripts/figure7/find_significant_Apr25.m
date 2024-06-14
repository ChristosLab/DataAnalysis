%find significant cells 
load('index_touse.mat');
load('all_activepassive_info.mat');
load('all_activepassive_data.mat');
all_location_names=unique(all_info(:,6));
group=all_info(:,3);  % pre post or error trials
%IndexC = strfind(group,'PASSIVE');
IndexC = strfind(group,'ACTIVE');
find_index = find(not(cellfun('isempty',IndexC)));
find_index=find_index(index_touse);
good_data=zeros(1,length(find_index)); %vectore indicate if cell was used

for i=1:length(find_index)
    
    cell_data=data_cell(find_index(i),:);
    cell_info=all_info(find_index(i),:);
    cell_location=cell_info(6);
    cell_location_num=strcmp(cell_location,all_location_names);
    vector_cell_loc(i)=find(cell_location_num==1); %1,Anterior-Dorsal 2, Anterior-Ventral 3, Mid-dorsal 4,Posterior-Dorsal 5,Posterior-Ventral
    
    all_fix=[];
    all_cue=[];
    all_cuedelay=[];
    all_sample=[];
    all_sampledelay=[];
    cue_loc=[];
    sample_loc=[];
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
        %two-way anova 
        if length(cueclass_ismatch)==length(cueclass_fix_count)
        %found some cells have ifmatch data but no spiking data
        good_data(i)=1;
        ifmatch=[ifmatch,cueclass_ismatch]; 
        temp_cue_loc=j*ones(1,length(cueclass_fix_count));
        cue_loc=[cue_loc,temp_cue_loc]; 
        temp_ismatch=cueclass_ismatch;
        sample_loc=[sample_loc,temp_cue_loc.*(temp_ismatch)+sample_loc_lookup(temp_cue_loc).*(temp_ismatch*(-1)+1)]; 
        all_fix=[all_fix,cueclass_fix_count];
        all_cue=[all_cue,cueclass_cue_rate];
        all_cuedelay=[all_cuedelay,cueclass_cuedelay_rate];
        all_sample=[all_sample,cueclass_sample_rate];
        all_sampledelay=[all_sampledelay,cueclass_sampledelay_rate];          
        end
    end
        cell_numtrials(i)=sum(cue_class_trials);

    % two-way anova for interaction, it is redundant with both one way and
    % twoway just to cross check results
    if (length(ifmatch)>0 && isnan(ifmatch(1))~=1)  %why some cell have ifmatch column NaN?
    [interaction_fix_p(i,:),tbl_fix]=anovan(all_fix,{cue_loc ifmatch},'model','interaction','display','off');
    f_fix(i,:)=cell2mat(tbl_fix(2:4,6))';
    [interaction_cue_p(i,:),tbl_cue]=anovan(all_cue,{cue_loc ifmatch},'model','interaction','display','off');
    f_cue(i,:)=cell2mat(tbl_cue(2:4,6))';
    [interaction_cuedelay_p(i,:),tbl_cuedelay]=anovan(all_cuedelay,{cue_loc ifmatch},'model','interaction','display','off');
    f_cuedelay(i,:)=cell2mat(tbl_cuedelay(2:4,6))';
    [interaction_sample_p(i,:),tbl_sample]=anovan(all_sample,{sample_loc ifmatch},'model','interaction','display','off');
    f_sample(i,:)=cell2mat(tbl_sample(2:4,6))';
    [interaction_sampledelay_p(i,:),tbl_sampledelay]=anovan(all_sampledelay,{sample_loc ifmatch},'model','interaction','display','off');
    f_sampledelay(i,:)=cell2mat(tbl_sampledelay(2:4,6))';
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
    find_sigfix_main2=intersect(find(interaction_fix_p(:,2)<0.05),data_inclusion_index);
    find_sigcue_main2=intersect(find(interaction_cue_p(:,2)<0.05),data_inclusion_index);
    find_sigcuedelay_main2=intersect(find(interaction_cuedelay_p(:,2)<0.05),data_inclusion_index);
    find_sigsample_main2=intersect(find(interaction_sample_p(:,2)<0.05),data_inclusion_index);
    find_sigsampledelay_main2=intersect(find(interaction_sampledelay_p(:,2)<0.05),data_inclusion_index);
    find_sigfix_interaction=intersect(find(interaction_fix_p(:,3)<0.05),data_inclusion_index);
    find_sigcue_interaction=intersect(find(interaction_cue_p(:,3)<0.05),data_inclusion_index);
    find_sigcuedelay_interaction=intersect(find(interaction_cuedelay_p(:,3)<0.05),data_inclusion_index);
    find_sigsample_interaction=intersect(find(interaction_sample_p(:,3)<0.05),data_inclusion_index);
    find_sigsampledelay_interaction=intersect(find(interaction_sampledelay_p(:,3)<0.05),data_inclusion_index);
    
    find_sigsample_CS1=setdiff(setdiff(find_sigsample_main1,find_sigsample_main2),find_sigsample_interaction);
    find_sigsample_CS2=setdiff(setdiff(find_sigsample_main2,find_sigsample_main1),find_sigsample_interaction); 
    find_sigsample_informative=unique([find_sigsample_main1;find_sigsample_main2;find_sigsample_interaction]);
    find_sigsample_LMS=setdiff(intersect(find_sigsample_main1,find_sigsample_main2),find_sigsample_interaction);
    find_sigsample_NMS=find_sigsample_interaction;
    find_sigsample_NMS1=intersect(find_sigsample_interaction,intersect(find_sigsample_main1,find_sigsample_main2));
    find_sigsample_NMS2=intersect(find_sigsample_interaction,setdiff(find_sigsample_main1,find_sigsample_main2));
    find_sigsample_NMS3=intersect(find_sigsample_interaction,setdiff(find_sigsample_main2,find_sigsample_main1));
    find_sigsample_NMS4=setdiff(find_sigsample_interaction,unique([find_sigsample_main1;find_sigsample_main2]));
    find_sigsample_CS=[find_sigsample_CS1;find_sigsample_CS2];
    
    find_sigsampledelay_CS1=setdiff(setdiff(find_sigsampledelay_main1,find_sigsampledelay_main2),find_sigsampledelay_interaction);
    find_sigsampledelay_CS2=setdiff(setdiff(find_sigsampledelay_main2,find_sigsampledelay_main1),find_sigsampledelay_interaction);
    find_sigsampledelay_informative=unique([find_sigsampledelay_main1;find_sigsampledelay_main2;find_sigsampledelay_interaction]);
    find_sigsampledelay_LMS=setdiff(intersect(find_sigsampledelay_main1,find_sigsampledelay_main2),find_sigsampledelay_interaction);
    find_sigsampledelay_NMS=find_sigsampledelay_interaction; 
    find_sigsampledelay_NMS1=intersect(find_sigsampledelay_interaction,intersect(find_sigsampledelay_main1,find_sigsampledelay_main2));
    find_sigsampledelay_NMS2=intersect(find_sigsampledelay_interaction,setdiff(find_sigsampledelay_main1,find_sigsampledelay_main2));
    find_sigsampledelay_NMS3=intersect(find_sigsampledelay_interaction,setdiff(find_sigsampledelay_main2,find_sigsampledelay_main1));
    find_sigsampledelay_NMS4=setdiff(find_sigsampledelay_interaction,unique([find_sigsampledelay_main1;find_sigsampledelay_main2]));
    find_sigsampledelay_CS=[find_sigsampledelay_CS1;find_sigsampledelay_CS2];
    
    prop_sigfix_main1=length(find_sigfix_main1)/length(data_inclusion_index);
    prop_sigcue_main1=length(find_sigcue_main1)/length(data_inclusion_index);
    prop_sigcuedelay_main1=length(find_sigcuedelay_main1)/length(data_inclusion_index);
    prop_sigsample_main1=length(find_sigsample_main1)/length(data_inclusion_index);
    prop_sigsampledelay_main1=length(find_sigsampledelay_main1)/length(data_inclusion_index);
    prop_sigfix_main2=length(find_sigfix_main2)/length(data_inclusion_index);
    prop_sigcue_main2=length(find_sigcue_main2)/length(data_inclusion_index);
    prop_sigcuedelay_main2=length(find_sigcuedelay_main2)/length(data_inclusion_index);
    prop_sigsample_main2=length(find_sigsample_main2)/length(data_inclusion_index);
    prop_sigsampledelay_main2=length(find_sigsampledelay_main2)/length(data_inclusion_index);
    prop_sigfix_interaction=length(find_sigfix_interaction)/length(data_inclusion_index);
    prop_sigcue_interaction=length(find_sigcue_interaction)/length(data_inclusion_index);
    prop_sigcuedelay_interaction=length(find_sigcuedelay_interaction)/length(data_inclusion_index);
    prop_sigsample_interaction=length(find_sigsample_interaction)/length(data_inclusion_index);
    prop_sigsampledelay_interaction=length(find_sigsampledelay_interaction)/length(data_inclusion_index);
    
    [find_informative,unique_ind]=unique([find_sigcue_main1;find_sigcuedelay_main1;find_sigsample_main1;find_sigsampledelay_main1;...  
    find_sigsample_main2;find_sigsampledelay_main2;...
    find_sigsample_interaction;find_sigsampledelay_interaction]);
    %proportion of informative cells in all recorded cells with good data,
    %informativ means main effect in any epoch of task
    proportion_informative=length(find_informative)/length(data_inclusion_index);
    
    figure;
    subplot(2,2,1);
    hold on;
    scatter(ones(1,length(data_inclusion_index)),f_fix(data_inclusion_index,1), 1);
    scatter(2*ones(1,length(data_inclusion_index)),f_cue(data_inclusion_index,1), 1);
    scatter(3*ones(1,length(data_inclusion_index)),f_cuedelay(data_inclusion_index,1), 1);
    scatter(4*ones(1,length(data_inclusion_index)),f_sample(data_inclusion_index,1), 1);
    scatter(5*ones(1,length(data_inclusion_index)),f_sampledelay(data_inclusion_index,1), 1);
    xticklabels({'fix','cue','cuedelay','sample','sampledelay'});
    ylim([0,100]);
    
    subplot(2,2,3);
    hold on;
    scatter(1,nanmean(f_fix(data_inclusion_index,1)), 50,'x');
    scatter(2,nanmean(f_cue(data_inclusion_index,1)), 50,'x');
    scatter(3,nanmean(f_cuedelay(data_inclusion_index,1)), 50,'x');
    scatter(4,nanmean(f_sample(data_inclusion_index,1)), 50,'x');
    scatter(5,nanmean(f_sampledelay(data_inclusion_index,1)), 50,'x');
    xticklabels({'fix','cue','cuedelay','sample','sampledelay'});  
    
    subplot(2,2,2);
    hold on;
    scatter(ones(1,length(data_inclusion_index)),f_fix(data_inclusion_index,2), 1);
    scatter(2*ones(1,length(data_inclusion_index)),f_cue(data_inclusion_index,2), 1);
    scatter(3*ones(1,length(data_inclusion_index)),f_cuedelay(data_inclusion_index,2), 1);
    scatter(4*ones(1,length(data_inclusion_index)),f_sample(data_inclusion_index,2), 1);
    scatter(5*ones(1,length(data_inclusion_index)),f_sampledelay(data_inclusion_index,2), 1);
    xticklabels({'fix','cue','cuedelay','sample','sampledelay'});
    ylim([0,60]);
    
    subplot(2,2,4);
    hold on;
    scatter(1,nanmean(f_fix(data_inclusion_index,2)), 80,'x');
    scatter(2,nanmean(f_cue(data_inclusion_index,2)), 80,'x');
    scatter(3,nanmean(f_cuedelay(data_inclusion_index,2)), 80,'x');
    scatter(4,nanmean(f_sample(data_inclusion_index,2)), 80,'x');
    scatter(5,nanmean(f_sampledelay(data_inclusion_index,2)), 80,'x');
    xticklabels({'fix','cue','cuedelay','sample','sampledelay'});  
    
    figure;
    main_effect_y=[prop_sigfix_main1,prop_sigcue_main1,prop_sigcuedelay_main1,prop_sigsample_main1,prop_sigsampledelay_main1;...
        prop_sigfix_main2,prop_sigcue_main2,prop_sigcuedelay_main2,prop_sigsample_main2,prop_sigsampledelay_main2;...
        prop_sigfix_interaction,prop_sigcue_interaction,prop_sigcuedelay_interaction,prop_sigsample_interaction,prop_sigsampledelay_interaction];
    bar(main_effect_y,'grouped');
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
    ylim([0,0.25]);
    text(1,0.2,strcat('n=',num2str(length(data_inclusion_index))));
    title('Active dataset');
    %title('Passive dataset');
    set(gca,'xticklabel',{'stim location','match/non-match','interaction'})
    legend('fix','cue','cuedelay','sample','sampledelay');
    figure;
    edges=[0:10:200];
    histogram(cell_numtrials,edges,'FaceAlpha',0.35);
    hold on;
    histogram(cell_numtrials(unique([find_sigcue_main1;find_sigsample_main1])),edges,'FaceAlpha',0.35);
    histogram(cell_numtrials(unique([find_sigsample_main2])),edges,'FaceAlpha',0.35);
 
    figure;
    labels={'CS','LMS','NMS'};
    subplot(1,2,1)
    mix_ratio_sample=[length(find_sigsample_CS)/length(find_sigsample_informative),...
        length(find_sigsample_LMS)/length(find_sigsample_informative),...
        length(find_sigsample_NMS)/length(find_sigsample_informative)];
    pie(mix_ratio_sample,labels);
    title(strcat('sample period n=',num2str(length(find_sigsample_informative))));
    subplot(1,2,2)
    mix_ratio_sampledelay=[length(find_sigsampledelay_CS)/length(find_sigsampledelay_informative),...
        length(find_sigsampledelay_LMS)/length(find_sigsampledelay_informative),...
        length(find_sigsampledelay_NMS)/length(find_sigsampledelay_informative)];
    pie(mix_ratio_sampledelay,labels);
    title(strcat('sampledelay period n=',num2str(length(find_sigsample_informative))));
    
    figure;
    labels={'CS','LMS','NMS','NS'};
    subplot(1,2,1)
    mix_ratio_sample=[length(find_sigsample_CS)/length(data_inclusion_index),...
        length(find_sigsample_LMS)/length(data_inclusion_index),...
        length(find_sigsample_NMS)/length(data_inclusion_index),...
        (length(data_inclusion_index)-length(find_sigsample_informative))/length(data_inclusion_index)];
    pie(mix_ratio_sample,labels);
    title(strcat('sample period n=',num2str(length(data_inclusion_index))));
    subplot(1,2,2)
    mix_ratio_sampledelay=[length(find_sigsampledelay_CS)/length(data_inclusion_index),...
        length(find_sigsampledelay_LMS)/length(data_inclusion_index),...
        length(find_sigsampledelay_NMS)/length(data_inclusion_index),...
        (length(data_inclusion_index)-length(find_sigsampledelay_informative))/length(data_inclusion_index)];
    pie(mix_ratio_sampledelay,labels);
    title(strcat('sampledelay period n=',num2str(length(data_inclusion_index))));
    disp(length(find_sigsample_CS1)/length(find_sigsample_CS));
    disp(length(find_sigsample_CS2)/length(find_sigsample_CS));
    disp(length(find_sigsampledelay_CS1)/length(find_sigsampledelay_CS));
    disp(length(find_sigsampledelay_CS2)/length(find_sigsampledelay_CS));
    disp(length(find_sigsample_NMS1)/length(find_sigsample_NMS));
    disp(length(find_sigsample_NMS2)/length(find_sigsample_NMS));
    disp(length(find_sigsample_NMS3)/length(find_sigsample_NMS));
    disp(length(find_sigsample_NMS4)/length(find_sigsample_NMS));
    disp(length(find_sigsampledelay_NMS1)/length(find_sigsampledelay_NMS));
    disp(length(find_sigsampledelay_NMS2)/length(find_sigsampledelay_NMS));
    disp(length(find_sigsampledelay_NMS3)/length(find_sigsampledelay_NMS));
    disp(length(find_sigsampledelay_NMS4)/length(find_sigsampledelay_NMS));
    
    index96=find(cell_numtrials>95);
    category_info{1}=unique([find_sigsample_NMS1]);
    category_info{2}=unique([find_sigsample_NMS2]);
    category_info{3}=unique([find_sigsample_NMS3]);
    category_info{4}=unique([find_sigsample_NMS4]);
    category_info{5}=unique([find_sigsample_CS1]);
    category_info{6}=unique([find_sigsample_CS2]);
    category_info{7}=unique([find_sigsample_LMS]);
   % save passive_sample_f.mat f_sample
   % save find_index.mat find_index
   % save active_trial96.mat index96
   % save active_informative.mat find_informative
   % save active_info.mat category_info
    %{
%plot specific brain area
figure;
plot_area=5;
area_name={'Anterior-Dorsal' 'Anterior-Ventral' 'Mid-dorsal' 'Posterior-Dorsal' 'Posterior-Ventral'};
all_areacells=intersect(find(vector_cell_loc==plot_area),data_inclusion_index);
area_fixmain1=intersect(find_sigfix_main1,find(vector_cell_loc==plot_area));
area_cuemain1=intersect(find_sigcue_main1,find(vector_cell_loc==plot_area));
area_cuedelaymain1=intersect(find_sigcuedelay_main1,find(vector_cell_loc==plot_area));
area_samplemain1=intersect(find_sigsample_main1,find(vector_cell_loc==plot_area));
area_sampledelaymain1=intersect(find_sigsampledelay_main1,find(vector_cell_loc==plot_area));
area_fixmain2=intersect(find_sigfix_main2,find(vector_cell_loc==plot_area));
area_cuemain2=intersect(find_sigcue_main2,find(vector_cell_loc==plot_area));
area_cuedelaymain2=intersect(find_sigcuedelay_main2,find(vector_cell_loc==plot_area));
area_samplemain2=intersect(find_sigsample_main2,find(vector_cell_loc==plot_area));
area_sampledelaymain2=intersect(find_sigsampledelay_main2,find(vector_cell_loc==plot_area));
area_fixinteraction=intersect(find_sigfix_interaction,find(vector_cell_loc==plot_area));
area_cueinteraction=intersect(find_sigcue_interaction,find(vector_cell_loc==plot_area));
area_cuedelayinteraction=intersect(find_sigcuedelay_interaction,find(vector_cell_loc==plot_area));
area_sampleinteraction=intersect(find_sigsample_interaction,find(vector_cell_loc==plot_area));
area_sampledelayinteraction=intersect(find_sigsampledelay_interaction,find(vector_cell_loc==plot_area));

prop_area_fixmain1=length(area_fixmain1)/length(all_areacells);
prop_area_cuemain1=length(area_cuemain1)/length(all_areacells);
prop_area_cuedelaymain1=length(area_cuedelaymain1)/length(all_areacells);
prop_area_samplemain1=length(area_samplemain1)/length(all_areacells);
prop_area_sampledelaymain1=length(area_sampledelaymain1)/length(all_areacells);
prop_area_fixmain2=length(area_fixmain2)/length(all_areacells);
prop_area_cuemain2=length(area_cuemain2)/length(all_areacells);
prop_area_cuedelaymain2=length(area_cuedelaymain2)/length(all_areacells);
prop_area_samplemain2=length(area_samplemain2)/length(all_areacells);
prop_area_sampledelaymain2=length(area_sampledelaymain2)/length(all_areacells);
prop_area_fixinteraction=length(area_fixinteraction)/length(all_areacells);
prop_area_cueinteraction=length(area_cueinteraction)/length(all_areacells);
prop_area_cuedelayinteraction=length(area_cuedelayinteraction)/length(all_areacells);
prop_area_sampleinteraction=length(area_sampleinteraction)/length(all_areacells);
prop_area_sampledelayinteraction=length(area_sampledelayinteraction)/length(all_areacells);

main_effect_area=[prop_area_fixmain1,prop_area_cuemain1,prop_area_cuedelaymain1,prop_area_samplemain1,prop_area_sampledelaymain1;...
        prop_area_fixmain2,prop_area_cuemain2,prop_area_cuedelaymain2,prop_area_samplemain2,prop_area_sampledelaymain2;...
         prop_area_fixinteraction,prop_area_cueinteraction,prop_area_cuedelayinteraction,prop_area_sampleinteraction,prop_area_sampledelayinteraction];
bar(main_effect_area,'grouped');
hold on;
ylim([0,0.35]);
xlimit = get(gca,'xlim');
plot(xlimit,[0.05 0.05]);
title(strcat(area_name{plot_area},' active dataset n=',num2str(length(all_areacells))));
%title(strcat(area_name{plot_area},' passive dataset n=',num2str(length(all_areacells))));

set(gca,'xticklabel',{'stim location','match/non-match','interaction'})
    legend('fix','cue','cuedelay','sample','sampledelay');
find_area_informative=intersect(find_informative,all_areacells); 
temp_index=find_index(find_area_informative);
%save post_dorsal_informative_post.mat temp_index
area_info{1}=all_areacells;
area_info{2}=find_area_informative;
area_info{3}=area_cuemain1;
area_info{4}=area_cuedelaymain1;
area_info{5}=area_samplemain1;
area_info{6}=area_sampledelaymain1;
area_info{7}=area_samplemain2;
area_info{8}=area_sampledelaymain2;
area_info{9}=find_index;
%}