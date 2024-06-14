%loop through resample from sample peudopopulation, different cell numbers
%sample with different size from randomly reassigned pool, result in the
%form of a curve, this is the number reported in the paper
load('post_sigsample_match.mat');
%load('pre_sigsampledelay_match.mat');
cell_select=temp;%informative cell index in all spatial cells

normalize=0;
load('all_spatial_info.mat');
load('all_spatial_data.mat');
%all_location_names=unique(all_spatial_info(:,5));
group=all_spatial_info(:,4);  
%IndexC = strfind(group,'PRE');
%IndexC = strfind(group,'POST');
%find_index = find(not(cellfun('isempty',IndexC)));
find_index=[1:length(group)];
good_data=zeros(1,length(find_index)); %vectore indicate if cell was used

for i=1:length(find_index)
    cell_data=all_spatial_data(find_index(i),:);
   % cell_info=all_spatial_info(find_index(i),:);
  %  cell_location=cell_info(5);
  %  cell_location_num=strcmp(cell_location,all_location_names);
  %  vector_cell_loc(i)=find(cell_location_num==1); %1,Anterior-Dorsal 2, Anterior-Ventral 3, Mid-dorsal 4,Posterior-Dorsal 5,Posterior-Ventral
    
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


edge1=2;
edge2=2.5;
cell_id=[];
cell_cueloc=[];
cell_ifmatch=[];
cell_bincount=[];
trial_count=0;

for i=1:length(data_inclusion_index) %loop cell

    temp_cell_data=all_spatial_data(data_inclusion_index(i),:);    
    for j=1:8  %loop cuecalss
        temp_cueclass_data=temp_cell_data{j};
        cell_id=[cell_id,i*ones(1,length(temp_cueclass_data))];
        cell_cueloc=[cell_cueloc,j*ones(1,length(temp_cueclass_data))];
        cell_ifmatch=[cell_ifmatch,[temp_cueclass_data.IsMatch]];
        cueclass_ontime=[temp_cueclass_data.Cue_onT];
        cueclass_spiketimes={temp_cueclass_data.TS};
        for p=1:length(temp_cueclass_data) %loop through trial
            trial_count=trial_count+1;
            temp_TS=cueclass_spiketimes{p}-cueclass_ontime(p);
            for q=1:length(edge1)  %loop through timebin            
            cell_bincount(trial_count,q)=length(find(temp_TS>edge1(q) & temp_TS<edge2(q)));
            end
        end
    end
end

%build sample pseudopopulation
for i=1:length(data_inclusion_index) %loop through cell
    disp(i);
    temp_cellindex=find(cell_id==i);
    ind_cell_bincount=cell_bincount(temp_cellindex,:);
    ind_cell_cueloc=cell_cueloc(temp_cellindex);
    ind_cell_ifmatch=cell_ifmatch(temp_cellindex);
    for j=1:8
        temp_index=find(ind_cell_ifmatch==1 & ind_cell_cueloc==j);
        rand_select=randperm(length(temp_index));
        select_index=rand_select(1:6);
     %   population_response(i,1+(j-1)*6:6+(j-1)*6,:)=ind_cell_bincount(temp_index(select_index),:);
        population_response(i,j,1:6)=ind_cell_bincount(temp_index(select_index),:);
    end
    for j=1:8
        temp_index=find(ind_cell_ifmatch==0 & ind_cell_cueloc==j);
        rand_select=randperm(length(temp_index));
        select_index=rand_select(1:6);
      %  population_response(i,48+1+(j-1)*6:48+6+(j-1)*6,:)=ind_cell_bincount(temp_index(select_index),:);
        population_response(i,8+j,1:6)=ind_cell_bincount(temp_index(select_index),:);
    end
end

if normalize==1
   for i=1:length(data_inclusion_index)
       temp_max=max(max(population_response(i,:,:)));
       if temp_max~=0
       population_response(i,:,:)=population_response(i,:,:)/temp_max;
       end
   end
end
%num_cell_vector=[100,200,300,500,1000,2000,3000,5000];
num_cell_vector=[200];
for n=1:length(num_cell_vector)
disp(n);
  %  if num_cell_vector(n)<size(population_response,1)
  %     all_winning_models=[];
  %     new_reponse=population_response;
  %     for i=1:10
  %     temp_index=randperm(size(new_response,1));
  %     test_response=new_response(temp_index(1:num_cell_vector(n)),:,:); %random select 100 cells
  %     wholebrain_all{1}=test_response;
  %     [subject_IDs, test_runs, winning_models, test_correlations]=functional_dimensionality_spike(wholebrain_all, ...
  %     ones(size(test_response,1),1),1,1);
  %     all_winning_models=[all_winning_models;winning_models];
  %     end
  %  else
  all_winning_models=[];
for j=1:50
       resample_times=ceil(num_cell_vector(n)/size(population_response,1));
       for i=1:resample_times*size(population_response,1)
       temp_pick=mod(i,size(population_response,1));
       if temp_pick==0
       temp_pick=size(population_response,1);
       end
       temp_response=squeeze(population_response(temp_pick,:,:));
       swap_match=randsample([0,1],1);
       spatial_resample=randperm(8);
       if swap_match==1
          reorder_temp_response=temp_response([spatial_resample+8,spatial_resample],:);
       else
          reorder_temp_response=temp_response([spatial_resample,spatial_resample+8],:);
       end
       new_response(i,:,:)=reorder_temp_response;
       end
       
       temp_index=randperm(size(new_response,1));
       test_response=new_response(temp_index(1:num_cell_vector(n)),:,:); %random select 100 cells
       wholebrain_all{1}=test_response;
       [subject_IDs, test_runs, winning_models, test_correlations]=functional_dimensionality_spike(wholebrain_all, ...
       ones(size(test_response,1),1),1,1);
      % all_winning_models=[all_winning_models;winning_models];
      all_winning_models=[all_winning_models,mean(winning_models)];
  %  end
end
mean_winning_model(n) = mean(all_winning_models);
ci_winning_model(n)=1.96*std(all_winning_models)/sqrt(30);
end

plot(mean_winning_model);