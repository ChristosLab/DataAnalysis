%loop through resample from different peudopopulation, different cell numbers
%sample from the same big pool

%load('post_sigstim_epoch.mat');
load('pre_sigdelay_epoch.mat'); %need to also change edge1 and edge2
cell_select=temp;%informative cell index in all feature cells


load('all_feature_info.mat');
load('all_feature_data.mat');
group=all_feature_info(:,4);  
%IndexC = strfind(group,'PRE');
%IndexC = strfind(group,'POST');
%find_index = find(not(cellfun('isempty',IndexC)));
find_index=[1:length(group)];
good_data=zeros(1,length(find_index)); %vectore indicate if cell was used

for i=1:length(find_index)
    cell_data=all_feature_data(find_index(i),:);
    
    for j=1:8 
        cueclass_data=cell_data{j};
        cue_class_trials(j)=length(cueclass_data);
        if cue_class_trials(j)>0
          if length(fieldnames(cueclass_data))>10
         %   cueclass_fix_count=[cueclass_data.fix];
            cueclass_ismatch=[cueclass_data.IsMatch];
          else
   %         cueclass_fix_count=[];
            cueclass_ismatch=[];  
          end
        cueclass_matchcount(1,j)=sum(cueclass_ismatch==1);
        cueclass_matchcount(2,j)=sum(cueclass_ismatch==0);
        cueclass_matchcount(3,j)=cueclass_matchcount(1,j)+cueclass_matchcount(2,j);
        else
        cueclass_matchcount(1,j)=0;
        cueclass_matchcount(2,j)=0;
        cueclass_matchcount(3,j)=0;
   %     cueclass_fix_count=[];
        cueclass_ismatch=[];
        end
        
    end
    cell_num_trial=sum(cue_class_trials);
    if cell_num_trial==sum(cueclass_matchcount(3,:))
       good_data(i)=1;
    end
    cueclass_quatify=(cueclass_matchcount(1,:)>5).*(cueclass_matchcount(2,:)>5);
    cell_quantify(i)=sum(cueclass_quatify);

end
data_inclusion_index=find_index(intersect(find(cell_quantify>7),find(good_data>0)));
data_inclusion_index=intersect(data_inclusion_index,cell_select);
%for predata set
%data_inclusion_index([])=[];


edge1=0.5;
edge2=2;
edge3=2.5;
edge4=4;
cell_id=[];
cell_stimloc=[];
cell_ifmatch=[];
cell_bincount=[];
cell_task_epoch=[];
trial_count=0;
sample_loc_lookup=[5 6 7 8 1 2 3 4]; %does not really match since only match trials are selected
for i=1:length(data_inclusion_index) %loop cell

    temp_cell_data=all_feature_data(data_inclusion_index(i),:);    
    for j=1:8  %loop cuecalss
        temp_cueclass_data=temp_cell_data{j};
        cell_id=[cell_id,i*ones(1,2*length(temp_cueclass_data))];
        cueclass_cueloc=j*ones(1,length(temp_cueclass_data));
        cell_ifmatch=[cell_ifmatch,[temp_cueclass_data.IsMatch],[temp_cueclass_data.IsMatch]];
        temp_ismatch=[temp_cueclass_data.IsMatch];
        cueclass_sampleloc=cueclass_cueloc.*(temp_ismatch)+sample_loc_lookup(cueclass_cueloc).*(temp_ismatch*(-1)+1);     
        cueclass_stimloc=[cueclass_cueloc,cueclass_sampleloc];
        cell_stimloc=[cell_stimloc,cueclass_stimloc];
        cell_task_epoch=[cell_task_epoch,[ones(1,length(temp_cueclass_data)),2*ones(1,length(temp_cueclass_data))]];

        cueclass_ontime=[temp_cueclass_data.Cue_onT];
        cueclass_spiketimes={temp_cueclass_data.TS};
        for p=1:length(temp_cueclass_data) %loop through trial
           % trial_count=trial_count+1;
            temp_TS=cueclass_spiketimes{p}-cueclass_ontime(p);
            cue_bincount=length(find(temp_TS>edge1 & temp_TS<edge2));
            cell_bincount=[cell_bincount,cue_bincount];
        end
        for p=1:length(temp_cueclass_data) %loop through trial
           % trial_count=trial_count+1;
            temp_TS=cueclass_spiketimes{p}-cueclass_ontime(p);
            sample_bincount=length(find(temp_TS>edge3 & temp_TS<edge4));
            cell_bincount=[cell_bincount,sample_bincount];
        end
    end
end


all_winning_models=[];

for p=1:50
%build sample pseudopopulation
for i=1:length(data_inclusion_index) %loop through cell
   % disp(i);
    temp_cellindex=find(cell_id==i);
    ind_cell_bincount=cell_bincount(temp_cellindex);
    ind_cell_stimloc=cell_stimloc(temp_cellindex);
    ind_cell_taskepoch=cell_task_epoch(temp_cellindex);
    ind_cell_ifmatch=cell_ifmatch(temp_cellindex);
    for j=1:8
        temp_index=find(ind_cell_ifmatch==1 & ind_cell_taskepoch==1 & ind_cell_stimloc==j);
        rand_select=randperm(length(temp_index));
        select_index=rand_select(1:6);
     %   population_response(i,1+(j-1)*6:6+(j-1)*6,:)=ind_cell_bincount(temp_index(select_index),:);
        population_response(i,j,1:6)=ind_cell_bincount(temp_index(select_index));
    end
    for j=1:8
        temp_index=find(ind_cell_ifmatch==1 & ind_cell_taskepoch==2 & ind_cell_stimloc==j);
        rand_select=randperm(length(temp_index));
        select_index=rand_select(1:6);
      %  population_response(i,48+1+(j-1)*6:48+6+(j-1)*6,:)=ind_cell_bincount(temp_index(select_index),:);
        population_response(i,8+j,1:6)=ind_cell_bincount(temp_index(select_index));
    end
end


temp_index=randperm(size(population_response,1));
test_response=population_response(temp_index(1:200),:,:); %random select 100 cells
wholebrain_all{1}=test_response;
[subject_IDs, test_runs, winning_models, test_correlations]=functional_dimensionality_spike(wholebrain_all, ...
    ones(size(test_response,1),1),1,1);
all_winning_models=[all_winning_models;mean(winning_models)];
end
mean_winning_model = mean(all_winning_models);
std_winning_model=1.96*std(all_winning_models)/sqrt(50);
