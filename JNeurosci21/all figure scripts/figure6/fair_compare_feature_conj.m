%this version uses 2 way anova for both conj and spatial dataset
load('all_feature_conj_info.mat');
load('new_conj_data.mat');
load('new_feature_data.mat');

unique_sets=unique(all_info(:,3));
feature10_group=[2,3,4,5,6,7];
featurex_group=[1,9,11,14,19,22,24];
featurexa_group=[8,10,12,15,16,18,20,23,25,26];
featurexb_group=[13,17,21,27,28];
for i=1:size(new_feature_cell,1)
    cell_set=all_info{i,3};
    cell_setnum=find(strcmp(unique_sets,cell_set));
    if ismember(cell_setnum,feature10_group)
        cell_group_num(i)=4;
    elseif ismember(cell_setnum,featurex_group)
        cell_group_num(i)=1;
    elseif ismember(cell_setnum,featurexa_group)
        cell_group_num(i)=2;
    else
        cell_group_num(i)=3;
    end
end

IndexC = [1:size(all_info,1)];
find_index = IndexC;
good_data_conj=zeros(1,length(find_index)); %vectore indicate if cell was used
good_data=zeros(1,length(find_index)); %for feature dataset

grouplookup{1}=[2,3,4,5,6,7,8,1];
grouplookup{2}=[2,3,8,5,6,7,4,1];
grouplookup{3}=[5,3,8,2,6,7,4,1];
grouplookup{4}=[2,3,4,5,6,7,8,1,10,9];
feature_dic=['C','D','H','N','P','S','T','Y'];

class_cueloc=[ones(1,8),2*ones(1,8),ones(1,8),2*ones(1,8)];
class_cuefeature=[ones(1,16),2*ones(1,16)];
class_sampleloc=repmat([1,1,2,2],1,8);
class_samplefeature=repmat([1,1,1,1,2,2,2,2],1,4);
%spatial_lookup=[5,6,7,8,1,2,3,4,5];
session_name=all_info(:,2);
cell_num=all_info(:,1);

for i=1:length(find_index)
    allconj_fix=[];
    allconj_cue=[];
    allconj_cuedelay=[];
    allconj_sample=[];
    allconj_sampledelay=[];
    allconj_cue_feature=[];
    allconj_sample_feature=[];
    allconj_ifmatch=[];
    allconj_cue_loc=[];
    allconj_sample_loc=[];
    conj_info=all_info(find_index(i),:);
    conj_set=conj_info{5};
     shape_to_use(1)=conj_set(1);
     shape1_code=find(feature_dic==shape_to_use(1));
     shape_to_use(2)=conj_set(2);
     shape2_code=find(feature_dic==shape_to_use(2));
    
    conj_data=new_conj_cell(find_index(i),:);
    conj_info=all_info(find_index(i),:);

%...............process conjunction..........................    
    for j=1:8  
        conj_cueclass_data=conj_data{j};
        conj_cueclass_trials(j)=length(conj_cueclass_data);
        if conj_cueclass_trials(j)>0
            if length(fieldnames(conj_cueclass_data))>10
        conj_fix_count=[conj_cueclass_data.fix];
        conj_cue_rate=[conj_cueclass_data.cuerate]; 
        conj_cuedelay_rate=[conj_cueclass_data.cuedelay];
        conj_sample_rate=[conj_cueclass_data.samplerate];
        conj_sampledelay_rate=[conj_cueclass_data.sampledelay];
        conj_ismatch=[conj_cueclass_data.IsMatch];
        conj_trialclass=[conj_cueclass_data.trialclass];
            else
        conj_fix_count=[];
        conj_cue_rate=[];
        conj_cuedelay_rate=[];
        conj_sample_rate=[];
        conj_sampledelay_rate=[];
        conj_ismatch=[];
        conj_trialclass=[];
            end
        else
        conj_fix_count=[];
        conj_cue_rate=[];
        conj_cuedelay_rate=[];
        conj_sample_rate=[];
        conj_sampledelay_rate=[];
        conj_ismatch=[];
        conj_trialclass=[];
        end
        %two-way anova 
        if length(conj_ismatch)==length(conj_fix_count)
        %found some cells have ifmatch data but no spiking data
        good_data_conj(i)=1;
        allconj_ifmatch=[allconj_ifmatch,conj_ismatch]; 
        
        temp_cue_feature=class_cuefeature(conj_trialclass);
        allconj_cue_feature=[allconj_cue_feature,temp_cue_feature];
        
        temp_cue_loc=class_cueloc(conj_trialclass);
        allconj_cue_loc=[allconj_cue_loc,temp_cue_loc]; 
        
        temp_sample_feature=class_samplefeature(conj_trialclass);
        allconj_sample_feature=[allconj_sample_feature,temp_sample_feature];
        
        temp_sample_loc=class_sampleloc(conj_trialclass);
        allconj_sample_loc=[allconj_sample_loc,temp_sample_loc];
        
        allconj_fix=[allconj_fix,conj_fix_count];
        allconj_cue=[allconj_cue,conj_cue_rate];
        allconj_cuedelay=[allconj_cuedelay,conj_cuedelay_rate];
        allconj_sample=[allconj_sample,conj_sample_rate];
        allconj_sampledelay=[allconj_sampledelay,conj_sampledelay_rate];          
        end
    end
    conj_select=find(allconj_cue_loc==allconj_sample_loc & allconj_cue_loc==1);
  %  conj_select=[1:length(allconj_cue_loc)];
%...............process conjunction..........................        
    cell_data=new_feature_cell(find_index(i),:);
    cell_info=all_info(find_index(i),:); 
    all_fix=[];  %for feature
    all_cue=[];
    all_cuedelay=[];
    all_sample=[];
    all_sampledelay=[];
    cue_feature=[];
    sample_feature=[];
    ifmatch=[];
    cell_feature_lookup=grouplookup{cell_group_num(find_index(i))};      
        
%...............process feature.............................. 
  for j=1:length(find(~cellfun(@isempty,cell_data)))      
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
        %found some cells have ifmatch data but no spiking data
        good_data(i)=1;
        ifmatch=[ifmatch,cueclass_ismatch]; 
        temp_cue_feature=j*ones(1,length(cueclass_fix_count));
        cue_feature=[cue_feature,temp_cue_feature]; 
        temp_ismatch=cueclass_ismatch;
        sample_feature=[sample_feature,temp_cue_feature.*(temp_ismatch)+cell_feature_lookup(temp_cue_feature).*(temp_ismatch*(-1)+1)]; 
        all_fix=[all_fix,cueclass_fix_count];
        all_cue=[all_cue,cueclass_cue_rate];
        all_cuedelay=[all_cuedelay,cueclass_cuedelay_rate];
        all_sample=[all_sample,cueclass_sample_rate];
        all_sampledelay=[all_sampledelay,cueclass_sampledelay_rate];          
        end
  end 
  feature_cue_select=find(cue_feature==shape1_code | cue_feature==shape2_code);
  %feature_cue_select=[1:length(cue_feature)];
  feature_sample_select=find(sample_feature==shape1_code | sample_feature==shape2_code);
  %feature_sample_select=[1:length(sample_feature)];
  %...............process feature..............................        
        
        
    cellcue_numtrials(i)=length(feature_cue_select); 
  %  cellsample_numtrials(i)=length(feature_sample_select);
    conj_numtrials(i)=length(conj_select);
   
    if (length(conj_ismatch)>0 && isnan(conj_ismatch(1))~=1)  %why some cell have ifmatch column NaN?
    [interaction_fix_conjp(i,:),tbl_fix]=anovan(allconj_fix(conj_select),{allconj_cue_feature(conj_select) allconj_ifmatch(conj_select)},'model','interaction','display','off');
    f_fix_conj(i,:)=cell2mat(tbl_fix(2:4,6))';
    [interaction_cue_conjp(i,:),tbl_cue]=anovan(allconj_cue(conj_select),{allconj_cue_feature(conj_select) allconj_ifmatch(conj_select)},'model','interaction','display','off');
    f_cue_conj(i,:)=cell2mat(tbl_cue(2:4,6))';
    [interaction_cuedelay_conjp(i,:),tbl_cuedelay]=anovan(allconj_cuedelay(conj_select),{allconj_cue_feature(conj_select) allconj_ifmatch(conj_select)},'model','interaction','display','off');
    f_cuedelay_conj(i,:)=cell2mat(tbl_cuedelay(2:4,6))';
    [interaction_sample_conjp(i,:),tbl_sample]=anovan(allconj_sample(conj_select),{allconj_sample_feature(conj_select) allconj_ifmatch(conj_select)},'model','interaction','display','off');
    f_sample_conj(i,:)=cell2mat(tbl_sample(2:4,6))';
    [interaction_sampledelay_conjp(i,:),tbl_sampledelay]=anovan(allconj_sampledelay(conj_select),{allconj_sample_feature(conj_select) allconj_ifmatch(conj_select)},'model','interaction','display','off');
    f_sampledelay_conj(i,:)=cell2mat(tbl_sampledelay(2:4,6))';
    else
    good_data_conj(i)=0;
    interaction_fix_conjp(i,1:3)=NaN;
    interaction_cue_conjp(i,1:3)=NaN;
    interaction_cuedelay_conjp(i,1:3)=NaN;
    interaction_sample_conjp(i,1:3)=NaN;
    interaction_sampledelay_conjp(i,1:3)=NaN;
    end
    
    if (length(ifmatch)>0 && isnan(ifmatch(1))~=1)  %why some cell have ifmatch column NaN?
    [interaction_fix_p(i,:),tbl_fix]=anovan(all_fix(feature_cue_select),{cue_feature(feature_cue_select) ifmatch(feature_cue_select)},'model','interaction','display','off');
    f_fix(i,:)=cell2mat(tbl_fix(2:4,6))';
    [interaction_cue_p(i,:),tbl_cue]=anovan(all_cue(feature_cue_select),{cue_feature(feature_cue_select) ifmatch(feature_cue_select)},'model','interaction','display','off');
    f_cue(i,:)=cell2mat(tbl_cue(2:4,6))';
    [interaction_cuedelay_p(i,:),tbl_cuedelay]=anovan(all_cuedelay(feature_cue_select),{cue_feature(feature_cue_select) ifmatch(feature_cue_select)},'model','interaction','display','off');
    f_cuedelay(i,:)=cell2mat(tbl_cuedelay(2:4,6))';
    [interaction_sample_p(i,:),tbl_sample]=anovan(all_sample(feature_sample_select),{sample_feature(feature_sample_select) ifmatch(feature_sample_select)},'model','interaction','display','off');
    f_sample(i,:)=cell2mat(tbl_sample(2:4,6))';
    [interaction_sampledelay_p(i,:),tbl_sampledelay]=anovan(all_sampledelay(feature_sample_select),{sample_feature(feature_sample_select) ifmatch(feature_sample_select)},'model','interaction','display','off');
    f_sampledelay(i,:)=cell2mat(tbl_sampledelay(2:4,6))';
    else
    good_data(i)=0;
    interaction_fix_p(i,1:3)=NaN;
    interaction_cue_p(i,1:3)=NaN;
    interaction_cuedelay_p(i,1:3)=NaN;
    interaction_sample_p(i,1:3)=NaN;
    interaction_sampledelay_p(i,1:3)=NaN;
    end
          
end

%...............plot conjunction results..................... 
    conj_inclusion_index=intersect(find(cellcue_numtrials>23),intersect(find(conj_numtrials>23),find(good_data_conj==1)));
    find_sigfix_main1=intersect(find(interaction_fix_conjp(:,1)<0.05),conj_inclusion_index);
    find_sigcue_main1=intersect(find(interaction_cue_conjp(:,1)<0.05),conj_inclusion_index);
    find_sigcuedelay_main1=intersect(find(interaction_cuedelay_conjp(:,1)<0.05),conj_inclusion_index);
    find_sigsample_main1=intersect(find(interaction_sample_conjp(:,1)<0.05),conj_inclusion_index);
    find_sigsampledelay_main1=intersect(find(interaction_sampledelay_conjp(:,1)<0.05),conj_inclusion_index);
    find_sigfix_main2=intersect(find(interaction_fix_conjp(:,2)<0.05),conj_inclusion_index);
    find_sigcue_main2=intersect(find(interaction_cue_conjp(:,2)<0.05),conj_inclusion_index);
    find_sigcuedelay_main2=intersect(find(interaction_cuedelay_conjp(:,2)<0.05),conj_inclusion_index);
    find_sigsample_main2=intersect(find(interaction_sample_conjp(:,2)<0.05),conj_inclusion_index);
    find_sigsampledelay_main2=intersect(find(interaction_sampledelay_conjp(:,2)<0.05),conj_inclusion_index);

    find_sigfix_interaction=intersect(find(interaction_fix_conjp(:,3)<0.05),conj_inclusion_index);
    find_sigcue_interaction=intersect(find(interaction_cue_conjp(:,3)<0.05),conj_inclusion_index);
    find_sigcuedelay_interaction=intersect(find(interaction_cuedelay_conjp(:,3)<0.05),conj_inclusion_index);
    find_sigsample_interaction=intersect(find(interaction_sample_conjp(:,3)<0.05),conj_inclusion_index);
    find_sigsampledelay_interaction=intersect(find(interaction_sampledelay_conjp(:,3)<0.05),conj_inclusion_index);

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
    
    prop_sigfix_main1=length(find_sigfix_main1)/length(conj_inclusion_index);
    prop_sigcue_main1=length(find_sigcue_main1)/length(conj_inclusion_index);
    prop_sigcuedelay_main1=length(find_sigcuedelay_main1)/length(conj_inclusion_index);
    prop_sigsample_main1=length(find_sigsample_main1)/length(conj_inclusion_index);
    prop_sigsampledelay_main1=length(find_sigsampledelay_main1)/length(conj_inclusion_index);
    prop_sigfix_main2=length(find_sigfix_main2)/length(conj_inclusion_index);
    prop_sigcue_main2=length(find_sigcue_main2)/length(conj_inclusion_index);
    prop_sigcuedelay_main2=length(find_sigcuedelay_main2)/length(conj_inclusion_index);
    prop_sigsample_main2=length(find_sigsample_main2)/length(conj_inclusion_index);
    prop_sigsampledelay_main2=length(find_sigsampledelay_main2)/length(conj_inclusion_index);

    prop_sigfix_interaction=length(find_sigfix_interaction)/length(conj_inclusion_index);
    prop_sigcue_interaction=length(find_sigcue_interaction)/length(conj_inclusion_index);
    prop_sigcuedelay_interaction=length(find_sigcuedelay_interaction)/length(conj_inclusion_index);
    prop_sigsample_interaction=length(find_sigsample_interaction)/length(conj_inclusion_index);
    prop_sigsampledelay_interaction=length(find_sigsampledelay_interaction)/length(conj_inclusion_index);

    [find_informative,unique_ind]=unique([find_sigcue_main1;find_sigcuedelay_main1;find_sigsample_main1;find_sigsampledelay_main1;...  
       find_sigcue_main2;find_sigcuedelay_main2;find_sigsample_main2;find_sigsampledelay_main2;...
       find_sigsample_interaction;find_sigsampledelay_interaction]);
   
    figure;
    main_effect_y=[
        prop_sigfix_main1,prop_sigcue_main1,prop_sigcuedelay_main1,prop_sigsample_main1,prop_sigsampledelay_main1;...
        prop_sigfix_main2,prop_sigcue_main2,prop_sigcuedelay_main2,prop_sigsample_main2,prop_sigsampledelay_main2;...
        prop_sigfix_interaction,prop_sigcue_interaction,prop_sigcuedelay_interaction,prop_sigsample_interaction,prop_sigsampledelay_interaction];
    bar(main_effect_y,'grouped');
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
    ylim([0,0.25]);
    text(1,0.2,strcat('n=',num2str(length(conj_inclusion_index))));
    title('Conjunction dataset');
    set(gca,'xticklabel',{'stim feature','match/non-match','feature x match'})
    legend('fix','cue','cuedelay','sample','sampledelay');
%...................plot conjunction results................
    figure;
    labels={'CS','LMS','NMS','NS'};
    subplot(1,2,1)
    mix_ratio_sample=[length(find_sigsample_CS)/length(conj_inclusion_index),...
        length(find_sigsample_LMS)/length(conj_inclusion_index),...
        length(find_sigsample_NMS)/length(conj_inclusion_index),...
        (length(conj_inclusion_index)-length(find_sigsample_informative))/length(conj_inclusion_index)];
    pie(mix_ratio_sample,labels);
    title(strcat('sample period n=',num2str(length(conj_inclusion_index))));
    subplot(1,2,2)
    mix_ratio_sampledelay=[length(find_sigsampledelay_CS)/length(conj_inclusion_index),...
        length(find_sigsampledelay_LMS)/length(conj_inclusion_index),...
        length(find_sigsampledelay_NMS)/length(conj_inclusion_index),...
        (length(conj_inclusion_index)-length(find_sigsampledelay_informative))/length(conj_inclusion_index)];
    pie(mix_ratio_sampledelay,labels);
    title(strcat('sampledelay period n=',num2str(length(conj_inclusion_index))));

     
    category_sample_info{1}=unique([find_sigsample_NMS1]);
    category_sample_info{2}=unique([find_sigsample_NMS2]);
    category_sample_info{3}=unique([find_sigsample_NMS3]);
    category_sample_info{4}=unique([find_sigsample_NMS4]);
    category_sample_info{5}=unique([find_sigsample_CS1]);
    category_sample_info{6}=unique([find_sigsample_CS2]);
    category_sample_info{7}=unique([find_sigsample_LMS]);
    category_sampledelay_info{1}=unique([find_sigsampledelay_NMS1]);
    category_sampledelay_info{2}=unique([find_sigsampledelay_NMS2]);
    category_sampledelay_info{3}=unique([find_sigsampledelay_NMS3]);
    category_sampledelay_info{4}=unique([find_sigsampledelay_NMS4]);
    category_sampledelay_info{5}=unique([find_sigsampledelay_CS1]);
    category_sampledelay_info{6}=unique([find_sigsampledelay_CS2]);
    category_sampledelay_info{7}=unique([find_sigsampledelay_LMS]);
    save find_index.mat find_index
    save conj_inclusion_index.mat conj_inclusion_index
    save conj_sample_f.mat f_sample_conj
    save conj_sampledelay_f.mat f_sampledelay_conj
    save conj_sample_info.mat category_sample_info
    save conj_sampledelay_info.mat category_sampledelay_info
%...................plot feature results....................
    data_inclusion_index=intersect(find(conj_numtrials>23),intersect(find(cellcue_numtrials>23),find(good_data==1)));
    
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
    title('feature dataset');
    set(gca,'xticklabel',{'stim feature','match/non-match','interaction'})
    legend('fix','cue','cuedelay','sample','sampledelay');
%...................plot spatial results....................

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
    
        category_sample_info{1}=unique([find_sigsample_NMS1]);
    category_sample_info{2}=unique([find_sigsample_NMS2]);
    category_sample_info{3}=unique([find_sigsample_NMS3]);
    category_sample_info{4}=unique([find_sigsample_NMS4]);
    category_sample_info{5}=unique([find_sigsample_CS1]);
    category_sample_info{6}=unique([find_sigsample_CS2]);
    category_sample_info{7}=unique([find_sigsample_LMS]);
    category_sampledelay_info{1}=unique([find_sigsampledelay_NMS1]);
    category_sampledelay_info{2}=unique([find_sigsampledelay_NMS2]);
    category_sampledelay_info{3}=unique([find_sigsampledelay_NMS3]);
    category_sampledelay_info{4}=unique([find_sigsampledelay_NMS4]);
    category_sampledelay_info{5}=unique([find_sigsampledelay_CS1]);
    category_sampledelay_info{6}=unique([find_sigsampledelay_CS2]);
    category_sampledelay_info{7}=unique([find_sigsampledelay_LMS]);
    save feature_inclusion_index.mat data_inclusion_index
    save feature_sample_f.mat f_sample
    save feature_sampledelay_f.mat f_sampledelay
    save feature_sample_info.mat category_sample_info
    save feature_sampledelay_info.mat category_sampledelay_info