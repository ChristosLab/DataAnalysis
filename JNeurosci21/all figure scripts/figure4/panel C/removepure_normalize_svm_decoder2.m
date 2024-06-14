function [x_time,time_perform,population_response]=removepure_normalize_svm_decoder(win_size,step_size,temp_index, mode,remove_mode,all_spatial_data)
% win_size and step_size to compute slide windown spike count,temp_index
% select cells to use, mode1 decode cue mode2 decode sample, mode3 decode
%remove_mode=1 remove matching pure, 2 remove location pure
% match/nonmatch
%this version convert 16 binary converter to try classification just by
%using output from fitcecoc (only one class code equles 1, all other 0).


%step_size=0.1;
%win_size=0.4;
%mode=3;
%remove_mode=1;
%temp_index=[2288,2276,1983,2001,2010];

if remove_mode==1
    condition_group=[1:8;9:16];
else
    condition_group=[1,9;2,10;3,11;4,12;5,13;6,14;7,15;8,16];
end
%load('all_spatial_data.mat');
%load('mid_dorsal_informative_pre.mat');
informative_select=temp_index;
normalize=1;  %if normalize the data

edge1=[-1:step_size:5-win_size];   %6 second of time 
edge2=[-1+win_size:step_size:5];
cell_id=[];
cell_cueloc=[];
cell_ifmatch=[];
cell_bincount=[];
store_num_trials=[];
trial_count=0;
%data extraction
for i=1:length(informative_select) %loop cell

    temp_cell_data=all_spatial_data(informative_select(i),:);    
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
match_label=[ones(1,48),zeros(1,48)];
cue_label=[];
sample_label=[];
condition_label=[];
cue_convert=[1:8,1:8];
sample_convert=[1:8,1:8];
match_convert=[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0];
lookup_table=[5,6,7,8,1,2,3,4];

for i=1:8
cue_label=[cue_label,i*ones(1,6)];
sample_label=[sample_label,i*ones(1,6)];
condition_label=[condition_label,i*ones(1,6)];
end
for i=1:8
cue_label=[cue_label,i*ones(1,6)];
sample_label=[sample_label,lookup_table(i)*ones(1,6)];
condition_label=[condition_label,(8+i)*ones(1,6)];
end

%build pseudopopulation
for i=1:length(informative_select) %loop through cell
    temp_cellindex=find(cell_id==i);
    ind_cell_bincount=cell_bincount(temp_cellindex,:);
    ind_cell_cueloc=cell_cueloc(temp_cellindex);
    ind_cell_ifmatch=cell_ifmatch(temp_cellindex);
    for j=1:8
        temp_index=find(ind_cell_ifmatch==1 & ind_cell_cueloc==j);
        rand_select=randperm(length(temp_index));
        select_index=rand_select(1:6);
        population_response(i,1+(j-1)*6:6+(j-1)*6,:)=ind_cell_bincount(temp_index(select_index),:);
    end
    for j=1:8
        temp_index=find(ind_cell_ifmatch==0 & ind_cell_cueloc==j);
        rand_select=randperm(length(temp_index));
        select_index=rand_select(1:6);
        population_response(i,48+1+(j-1)*6:48+6+(j-1)*6,:)=ind_cell_bincount(temp_index(select_index),:);
    end
end

if normalize==1
for x=1:size(population_response,1)
    temp_data=squeeze(population_response(x,:,:));
    cell_mean=mean(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    cell_std=std(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    population_response(x,:,:)=(population_response(x,:,:)-cell_mean)./cell_std;
end
end

%remove pure selectivity
for i=1:size(population_response,1)
  for t=1:size(population_response,3)
    cell_timepoint_response=population_response(i,:,t);
    cell_response=mat2cell(reshape(cell_timepoint_response,6,16)',ones(1,16),6);
    removepure_response(i,:,t)=cell2mat(remove_pure(cell_response',condition_group));
  end
end
%population_response=removepure_response;
%{
num_cell_target=[450];
resample_times=ceil(num_cell_target/size(population_response,1));
for i=1:resample_times*size(population_response,1)
   temp_pick=mod(i,size(population_response,1));
   if temp_pick==0
       temp_pick=size(population_response,1);
   end
   temp_response=squeeze(population_response(temp_pick,:,:));
   swap_match=randsample([0,1],1);
   temp_resample=randperm(8);
   spatial_resample=[(temp_resample-1)*6+1;(temp_resample-1)*6+2;(temp_resample-1)*6+3;(temp_resample-1)*6+4;(temp_resample-1)*6+5;(temp_resample-1)*6+6];
   spatial_resample=reshape(spatial_resample,1,48);
   if swap_match==1
      reorder_temp_response=temp_response([spatial_resample+48,spatial_resample],:);
   else
      reorder_temp_response=temp_response([spatial_resample,spatial_resample+48],:);
   end
   new_response(i,:,:)=reorder_temp_response;
end
population_response=new_response(1:num_cell_target,:,:);
%}
%t=templateSVM('KernelFunction','Linear');
t=templateSVM('KernelFunction','polynomial','polynomialOrder',2,'BoxConstraint',1);

if mode==1
   true_decode_label=cue_label;
elseif mode==2
   true_decode_label=sample_label;
else
   true_decode_label=match_label;
end

decode_label=condition_label;
  %  decode_label=decode_label(randperm(length(decode_label)));
for time_point=1:length(edge1)
    disp(time_point);

OBJ=fitcecoc(population_response(:,:,time_point)',decode_label,'Coding','onevsall','Learners',t);
cvecoc = crossval(OBJ,'KFold',5);
%cvecoc = crossval(OBJ,'leaveout','on');
Yhat = kfoldPredict(cvecoc);
if mode==1
   true_Yhat=cue_convert(Yhat);
elseif mode==2
   true_Yhat=sample_convert(Yhat);
else    
   true_Yhat=match_convert(Yhat);
end
confuse_mat=confusionmat(true_decode_label,true_Yhat);
num_correct=diag(confuse_mat);
correct_pertrue_class=num_correct./sum(confuse_mat,2);
mean_correct_linearSVM=mean(correct_pertrue_class);
time_perform(time_point)=mean_correct_linearSVM;
crv_cell{time_point}=cvecoc;
end
x_time=linspace(-1,5,length(edge1))+0.5*win_size;
end
