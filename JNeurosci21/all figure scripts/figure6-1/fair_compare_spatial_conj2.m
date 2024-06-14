%this version uses 2 way anova for both conj and spatial dataset
load('all_feature_conj_info.mat');
load('new_conj_data.mat');
load('all_spatial_info.mat');
load('all_spatial_data.mat');
all_location_names=unique(all_info(:,7));
mode=1; %mode1 use class 1 and 3 in conj data, mode2 use classe 6 and 8 in conj data
if mode==1
class_touse=[1,3];
else
class_touse=[6,8];
end
IndexC = [1:size(all_info,1)];
find_index = IndexC;
good_data_conj=zeros(1,length(find_index)); %vectore indicate if cell was used
good_data=zeros(1,length(find_index));
class_cueloc=[ones(1,8),2*ones(1,8),ones(1,8),2*ones(1,8)];
class_cuefeature=[ones(1,16),2*ones(1,16)];
class_sampleloc=repmat([1,1,2,2],1,8);
class_samplefeature=repmat([1,1,1,1,2,2,2,2],1,4);
spatial_lookup=[5,6,7,8,1,2,3,4,5];
session_name=all_info(:,2);
cell_num=all_info(:,1);
for i=1:length(session_name)
    temp_session_name=session_name{i};
    temp_session_name=temp_session_name(1:end-2);
    temp_cell_name=cell_num{i};
    spatial_search{i}=[temp_session_name,'_1_',temp_cell_name];
end
for i=1:length(all_spatial_info(:,1))  %get all spatial dataset names
    temp_name1=all_spatial_info{i,1};
    temp_name2=num2str(all_spatial_info{i,2});
    spatial_names{i}=[temp_name1(1:8),'_',temp_name2];
end
for i=1:length(spatial_search) %search if paired spatial file exist
    if length(find(ismember(spatial_names,spatial_search{i})==1))>0
    ind(i)=find(ismember(spatial_names,spatial_search{i})==1);
    else
    ind(i)=NaN;    
    end
end
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
    loc_to_use(1)=str2num(conj_set(3));
    
    if (isnan(ind(i))==0 & loc_to_use(1)~=9)%judge if conjunction data uses 9 as first location
    conj_data=new_conj_cell(find_index(i),:);
    conj_info=all_info(find_index(i),:);
    spatial_data=all_spatial_data(ind(i),:);
    spatial_info=all_spatial_info(ind(i),:);
    cell_location=conj_info(7);
    cell_location_num=strcmp(cell_location,all_location_names);
    vector_cell_loc(i)=find(cell_location_num==1); %1,Anterior-Dorsal 2, Anterior-Ventral 3, Mid-dorsal 4,Posterior-Dorsal 5,Posterior-Ventral   
    conj_set=conj_info{5};
    loc_to_use(1)=str2num(conj_set(3));
    loc_to_use(2)=spatial_lookup(loc_to_use(1));
    
%...............process conjunction..........................    
    for j=1:2  %use only 2 classes for conjunction
        conj_cueclass_data=conj_data{class_touse(j)};
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
%...............process conjunction..........................        
     
    all_fix=[];   %for spatial
    all_cue=[];
    all_cuedelay=[];
    all_sample=[];
    all_sampledelay=[];
    ifmatch=[];
    cue_loc=[];
    sample_loc=[];       
        
%...............process spatial.............................. 
    for j=1:2
        cueclass_data=spatial_data{loc_to_use(j)};
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
        temp_cue_loc=loc_to_use(j)*ones(1,length(cueclass_fix_count));
        cue_loc=[cue_loc,temp_cue_loc]; 
        temp_ismatch=cueclass_ismatch;
        sample_loc=[sample_loc,temp_cue_loc.*(temp_ismatch)+spatial_lookup(temp_cue_loc).*(temp_ismatch*(-1)+1)]; 
        all_fix=[all_fix,cueclass_fix_count];
        all_cue=[all_cue,cueclass_cue_rate];
        all_cuedelay=[all_cuedelay,cueclass_cuedelay_rate];
        all_sample=[all_sample,cueclass_sample_rate];
        all_sampledelay=[all_sampledelay,cueclass_sampledelay_rate];          
        end
    end
%...............process spatial..............................        
        
        
    cell_numtrials(i)=sum(cue_class_trials);   
    conj_numtrials(i)=sum(conj_cueclass_trials);
   
    if (length(conj_ismatch)>0 && isnan(conj_ismatch(1))~=1)  %why some cell have ifmatch column NaN?
    [interaction_fix_conjp(i,:),tbl_fix]=anovan(allconj_fix,{allconj_cue_loc allconj_ifmatch},'model','interaction','display','off');
    f_fix_conj(i,:)=cell2mat(tbl_fix(2:4,6))';
    [interaction_cue_conjp(i,:),tbl_cue]=anovan(allconj_cue,{allconj_cue_loc allconj_ifmatch},'model','interaction','display','off');
    f_cue_conj(i,:)=cell2mat(tbl_cue(2:4,6))';
    [interaction_cuedelay_conjp(i,:),tbl_cuedelay]=anovan(allconj_cuedelay,{allconj_cue_loc allconj_ifmatch},'model','interaction','display','off');
    f_cuedelay_conj(i,:)=cell2mat(tbl_cuedelay(2:4,6))';
    [interaction_sample_conjp(i,:),tbl_sample]=anovan(allconj_sample,{allconj_sample_loc allconj_ifmatch},'model','interaction','display','off');
    f_sample_conj(i,:)=cell2mat(tbl_sample(2:4,6))';
    [interaction_sampledelay_conjp(i,:),tbl_sampledelay]=anovan(allconj_sampledelay,{allconj_sample_loc allconj_ifmatch},'model','interaction','display','off');
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
    end
      
    end
end

%...............plot conjunction results..................... 
    conj_inclusion_index=intersect(find(cell_numtrials>23),intersect(find(conj_numtrials>23),find(good_data_conj==1)));
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
    ylim([0,0.4]);
    text(1,0.2,strcat('n=',num2str(length(conj_inclusion_index))));
    title('Conjunction dataset');
    set(gca,'xticklabel',{'stim location','match/non-match','loc x match'})
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
%...................plot spatial results....................
    data_inclusion_index=intersect(find(conj_numtrials>23),intersect(find(cell_numtrials>23),find(good_data==1)));
    
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
    ylim([0,0.4]);
    text(1,0.2,strcat('n=',num2str(length(data_inclusion_index))));
    title('spatial dataset');
    set(gca,'xticklabel',{'stim location','match/non-match','interaction'})
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
    save spatial_inclusion_index.mat data_inclusion_index
    save spatial_sample_f.mat f_sample
    save spatial_sampledelay_f.mat f_sampledelay
    save spatial_sample_info.mat category_sample_info
    save spatial_sampledelay_info.mat category_sampledelay_info