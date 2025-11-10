function [x_time,time_perform,crv_cell,population_response]=linear_decoder(win_size,step_size,temp_index, mode)
% win_size and step_size to compute slide windown spike count,temp_index
% select cells to use, mode1 decode cue mode2 decode sample, mode3 decode
% match/nonmatch
load('test_data.mat');
%load('mid_dorsal_informative_pre.mat');
informative_select=temp_index;
normalize=1;  %if normalize the data
pca_denoice=1; %if pca denoise the data
%step_size=0.1;
%win_size=0.2;
edge1=[-1.5:step_size:5-win_size];   %6 second of time 
edge2=[-1.5+win_size:step_size:5];

%edge1 = 0.75;
%edge2 = 1.75; 
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
lookup_table=[2,3,4,5,6,7,8,1];

for i=1:8
cue_label=[cue_label,i*ones(1,6)];
sample_label=[sample_label,i*ones(1,6)];
end
for i=1:8
cue_label=[cue_label,i*ones(1,6)];
sample_label=[sample_label,lookup_table(i)*ones(1,6)];
end

% 6 x8 = 48, 48 x2 = 96
%build pseudopopulation
for i=1:length(informative_select) %loop through cell
    temp_cellindex=find(cell_id==i);
    ind_cell_bincount=cell_bincount(temp_cellindex,:);
    ind_cell_cueloc=cell_cueloc(temp_cellindex);
    ind_cell_ifmatch=cell_ifmatch(temp_cellindex);
    if length(find(ind_cell_ifmatch==1))>47 & length(find(ind_cell_ifmatch==0))>47
    for j=1:8
        temp_index=find(ind_cell_ifmatch==1 & ind_cell_cueloc==j);
        rand_select=randperm(length(temp_index)); %only has length of 6, (12/2)
        if length(rand_select)<6
            disp('debug');
        end
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
end
%end

if normalize==1
    zero_count=0;
    zero_index=[];
for x=1:size(population_response)
    temp_data=squeeze(population_response(x,:,:));
    cell_mean=mean(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    cell_std=std(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    population_response(x,:,:)=(population_response(x,:,:)-cell_mean)./cell_std;
    if cell_std==0
    zero_count=zero_count+1;
    zero_index=[zero_index,x];
    end
end
    population_response(zero_index,:,:)=[];
end
%if pcd_denoise==1
%component_touse=3;
%[coeff,score,latent] = pca(population_response);
%end
t=templateLinear();
%svm for each time point
for i=1:length(edge1)
    %disp(i);
    if mode==1
        decode_label=cue_label;
    elseif mode==2
        decode_label=sample_label;
    else
        decode_label=match_label;
    end
    
    %OBJ=fitcecoc(population_response(:,:,i)',decode_label,'Learners',t);
    %OBJ=fitcecoc(population_response(:,:,i)',decode_label,"Learners",'Linear');
    OBJ=fitcecoc(population_response(:,:,i)',decode_label,'Learners','linear','CrossVal','on', 'KFold',10);
    cvecoc=OBJ;
    %OBJ=fitcdiscr(population_response(:,:,i)',decode_label);
    %OBJ=fitcdiscr(population_response(:,:,i)',decode_label, 'Coding','onevsall');
    % pop response 78x110x57
    % decode label 1X96
    
    %cvecoc = crossval(OBJ,'KFold',5);
   % cvecoc = crossval(OBJ,'leaveout','on');
    Yhat = kfoldPredict(cvecoc);
    confuse_mat=confusionmat(decode_label,Yhat);
    num_correct=diag(confuse_mat);
    correct_pertrue_class=num_correct./sum(confuse_mat,2);
    mean_correct_linearSVM=mean(correct_pertrue_class);
    time_perform(i)=mean_correct_linearSVM;
    crv_cell{i}=cvecoc;
end
x_time=linspace(-1.5,5,length(edge1))+0.5*win_size;
end