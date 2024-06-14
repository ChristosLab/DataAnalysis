function [time_perform,population_response]=numbi_normalize_svm_decoder(temp_index,all_spatial_data)
% win_size and step_size to compute slide windown spike count,temp_index
% select cells to use, mode1 decode cue mode2 decode sample, mode3 decode
%remove_mode=1 remove matching pure, 2 remove location pure
% this version use tru error correcting code

%step_size=0.1;
%win_size=0.4;

informative_select=temp_index;
normalize=0;  %if normalize the data

edge1=[2];   %6 second of time 
edge2=[2.5];
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
%{
% code to remove pure    
   condition_group=[1,9;2,10;3,11;4,12;5,13;6,14;7,15;8,16];
for i=1:size(population_response,1)
    cell_timepoint_response=population_response(i,:);
    cell_response=mat2cell(reshape(cell_timepoint_response,6,16)',ones(1,16),6);
    removepure_response(i,:)=cell2mat(remove_pure(cell_response',condition_group));
end
population_response=removepure_response;
%}
%{
%code to radomly reconstruct
       num_cell_vector=[200];
       resample_times=ceil(num_cell_vector/size(population_response,1));
       for i=1:resample_times*size(population_response,1)
       temp_pick=mod(i,size(population_response,1));
       if temp_pick==0
       temp_pick=size(population_response,1);
       end
       temp_response=squeeze(population_response(temp_pick,:));
       swap_match=randsample([0,1],1);
       spatial_resample=repmat(randperm(8)*6,6,1)+repmat([-5,-4,-3,-2,-1,0]',1,8);
       spatial_resample=reshape(spatial_resample,1,48);
       if swap_match==1
          reorder_temp_response=temp_response([spatial_resample+48,spatial_resample]);
       else
          reorder_temp_response=temp_response([spatial_resample,spatial_resample+48]);
       end
       new_response(i,:)=reorder_temp_response;
       end
population_response=new_response;
%}
t=templateSVM('KernelFunction','Linear');
%t=templateSVM('KernelFunction','polynomial','polynomialOrder',2,'BoxConstraint',1);

%all_labels=[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0];
all_labels=dec2bin(sum(nchoosek(2.^(0:16-1),8),2)) - '0';
hundred_count=0;
old_count=0;
for binary_classification=1:size(all_labels,1)  
decode_label=reshape(repmat(all_labels(binary_classification,:),6,1),1,96);
OBJ=fitcecoc(population_response(:,:)',decode_label,'Coding','onevsall','Learners',t);
cvecoc = crossval(OBJ,'KFold',5);
hundred_count=ceil(binary_classification/100);
if hundred_count>old_count
   disp(hundred_count);
   if hundred_count==30
       disp('test');
   end
end
old_count=hundred_count;
Yhat = kfoldPredict(cvecoc);
confuse_mat=confusionmat(decode_label,Yhat);
num_correct=diag(confuse_mat);
correct_pertrue_class=num_correct./sum(confuse_mat,2);
mean_correct_linearSVM(binary_classification)=mean(correct_pertrue_class);
end

time_perform=mean_correct_linearSVM;
%crv_cell{time_point}=cvecoc;

end
