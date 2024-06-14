function [x_time,time_perform,crv_cell,population_response]=SVM_decoder3(win_size,step_size,temp_index, mode,stimset)
% this version works with feature datasets
load('all_feature_data.mat');
normalize=1;  %if normalize the data
pca_denoice=1; %if pca denoise the data
grouplookup{1}=[2,3,4,5,6,7,8,1];
grouplookup{2}=[2,3,8,5,6,7,4,1];
grouplookup{3}=[5,3,8,2,6,7,4,1];
grouplookup{4}=[2,3,4,5,6,7,8,1,10,9];
lookup_table=grouplookup{stimset};

informative_select=temp_index;

edge1=[-1:step_size:5-win_size];   %6 second of time 
edge2=[-1+win_size:step_size:5];
cell_id=[];
cell_cueloc=[];
cell_ifmatch=[];
cell_bincount=[];
trial_count=0;
good_index=ones(1,length(informative_select));
%data extraction
for i=1:length(informative_select)
    temp_cell_data=all_feature_data(informative_select(i),:);   
    for j=1:8
        temp_cueclass_data=temp_cell_data{j};
        cueclass_numtrial(j)=length(temp_cueclass_data);
    end 
    if min(cueclass_numtrial)<12
        good_index(i)=0;
    end
end
new_select=informative_select(good_index==1);

for i=1:length(new_select) %loop cell
    temp_cell_data=all_feature_data(new_select(i),:);    
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
for p=1:8
  cue_label=[cue_label,p*ones(1,6)];
  sample_label=[sample_label,p*ones(1,6)];
end
for p=1:8
  cue_label=[cue_label,p*ones(1,6)];
  sample_label=[sample_label,lookup_table(p)*ones(1,6)];
end
    
%build pseudopopulation
for i=1:length(new_select) %loop through cell
   % disp(i);   
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
for x=1:size(population_response)
    temp_data=squeeze(population_response(x,:,:));
    cell_mean=mean(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    cell_std=std(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    population_response(x,:,:)=(population_response(x,:,:)-cell_mean)./cell_std;
end
end

t=templateSVM('KernelFunction','linear');
%svm for each time point
for i=1:length(edge1)
   % disp(i);
    if mode==1
        decode_label=cue_label;
    elseif mode==2
        decode_label=sample_label;
    else
        decode_label=match_label;
    end
    OBJ=fitcecoc(population_response(:,:,i)',decode_label,'Coding','onevsall','Learners',t);
    cvecoc = crossval(OBJ,'KFold',10);
   % cvecoc = crossval(OBJ,'leaveout','on');
    Yhat = kfoldPredict(cvecoc);
    confuse_mat=confusionmat(decode_label,Yhat);
    num_correct=diag(confuse_mat);
    correct_pertrue_class=num_correct./sum(confuse_mat,2);
    mean_correct_linearSVM=mean(correct_pertrue_class);
    time_perform(i)=mean_correct_linearSVM;
    crv_cell{i}=cvecoc;
end
x_time=linspace(-1,5,length(edge1))+0.5*win_size;
end
