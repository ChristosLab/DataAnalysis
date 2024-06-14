function [x_time,time_perform_correct,time_perform_err,time_perform_all]=spatial_saccade_decoder(win_size,step_size,correctcell_index,errcell_index,condition_index)
%this decoder decod both correct and err trials
load('all_spatial_data.mat');
load('all_spatialerr_data.mat');
load('spatial_sacdirect.mat');
%load('mid_dorsal_informative_pre.mat');
informative_select=correctcell_index;
err_informative_select=errcell_index;
normalize=0;  %if normalize the data
num_trial_percondition=2;
%step_size=0.1;
%win_size=0.2;
edge1=[-1:step_size:6-win_size];   %6 second of time 
edge2=[-1+win_size:step_size:6];
correct_cell_id=[];
correct_ifmatch=[];
correct_bincount=[];
correct_sacdirect=[];
err_cell_id=[];
err_ifmatch=[];
err_bincount=[];
err_sacdirect=[];
%store_num_trials=[];
correct_trial_count=0;
err_trial_count=0;
%data extraction
for i=1:length(informative_select) %loop cell
    cell_saccade=[all_trialdirect{informative_select(i)}]; %all saccade direction of the cell
    correct_cell_data=all_spatial_data(informative_select(i),:);  
    err_cell_data=all_spatialerr_data(err_informative_select(i),:);    
    for j=1:8  %loop cuecalss
        temp_correctdata=correct_cell_data{j};%get physioogy data
        correct_trialnum=[temp_correctdata.trialnum]; %correct trial index
        correct_class_sacdirect=cell_saccade(correct_trialnum); %correct trial sac direction
        correct_cell_id=[correct_cell_id,i*ones(1,length(temp_correctdata))]; %store task information
        correct_ifmatch=[correct_ifmatch,[temp_correctdata.IsMatch]];
        correct_sacdirect=[correct_sacdirect,correct_class_sacdirect];
        cueclass_ontime=[temp_correctdata.Cue_onT]; %store spike information
        cueclass_spiketimes={temp_correctdata.TS};
        for p=1:length(temp_correctdata) %loop through correct trials
            correct_trial_count=correct_trial_count+1;
            temp_TS=cueclass_spiketimes{p}-cueclass_ontime(p);
            for q=1:length(edge1)  %loop through timebin            
            correct_bincount(correct_trial_count,q)=length(find(temp_TS>edge1(q) & temp_TS<edge2(q)));
            end
        end
        
       
        temp_errdata=err_cell_data{j};
        if length(temp_errdata)>0
        if isfield(temp_errdata,'trialnum')
        err_trialnum=[temp_errdata.trialnum];
        else
        err_trialnum=[temp_errdata.Trial_Num];
        end
        err_class_sacdirect=cell_saccade(err_trialnum);
        err_cell_id=[err_cell_id,i*ones(1,length(temp_errdata))]; %store task information
        err_ifmatch=[err_ifmatch,[temp_errdata.IsMatch]];
        err_sacdirect=[err_sacdirect,err_class_sacdirect];       
        cueclass_ontime=[temp_errdata.Cue_onT]; %store spike information
        cueclass_spiketimes={temp_errdata.TS};
        for p=1:length(temp_errdata) %loop through correct trials
            err_trial_count=err_trial_count+1;
            temp_TS=cueclass_spiketimes{p}-cueclass_ontime(p);
            for q=1:length(edge1)  %loop through timebin            
            err_bincount(err_trial_count,q)=length(find(temp_TS>edge1(q) & temp_TS<edge2(q)));
            end
        end
        end
        
    end
end
%build peudopopulaiton for correct data
for i=1:length(informative_select) %loop through cell
    correct_cellindex=find(correct_cell_id==i); 
    ind_cell_bincount=correct_bincount(correct_cellindex,:);%construct corret cell population response
   % ind_cell_cueloc=correct_cueloc(correct_cellindex);
    ind_cell_ifmatch=correct_ifmatch(correct_cellindex);
    ind_cell_sacdirect=correct_sacdirect(correct_cellindex);
    for j=1:8
        temp_index=find(ind_cell_sacdirect==j & ind_cell_ifmatch==1);
        rand_select=randperm(length(temp_index));
        select_index=rand_select(1:6);
        correct_population_response(i,1+(j-1)*6:6+(j-1)*6,:)=ind_cell_bincount(temp_index(select_index),:);
        temp_index=find(ind_cell_sacdirect==j & ind_cell_ifmatch==0);
        rand_select=randperm(length(temp_index));
        select_index=rand_select(1:6);
        correct_population_response(i,49+(j-1)*6:54+(j-1)*6,:)=ind_cell_bincount(temp_index(select_index),:);
    end
    
end
%build pseudopopulation for err data
for i=1:length(informative_select) %loop through cell
    
    err_cellindex=find(err_cell_id==i); 
    ind_cell_bincount=err_bincount(err_cellindex,:);%construct err cell population response
    ind_cell_ifmatch=err_ifmatch(err_cellindex);
    ind_cell_sacdirect=err_sacdirect(err_cellindex);
    sac_err_label=[];
    sac_correct_label=[];
    for j=1:length(condition_index)
        current_condition=condition_index(j);
        current_match=(current_condition<9);
        current_sacdirect=mod(current_condition,8);
        if current_sacdirect==0
           current_sacdirect=8;
        end
        sac_err_label=[sac_err_label,current_sacdirect*ones(1,num_trial_percondition)];
        sac_correct_label=[sac_correct_label,current_sacdirect*ones(1,6)];
        temp_index=find(ind_cell_ifmatch==current_match & ind_cell_sacdirect==current_sacdirect);
        rand_select=randperm(length(temp_index));
        select_index=rand_select(1:num_trial_percondition);
        err_population_response(i,1+(j-1)*num_trial_percondition:j*num_trial_percondition,:)=ind_cell_bincount(temp_index(select_index),:);
    end
end
%end
if normalize==1
for x=1:size(correct_population_response)
    temp_data=squeeze(correct_population_response(x,:,:));
    cell_mean=mean(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    cell_std=std(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    correct_population_response(x,:,:)=(correct_population_response(x,:,:)-cell_mean)./cell_std;
end
for x=1:size(err_population_response)
    temp_data=squeeze(err_population_response(x,:,:));
    cell_mean=mean(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    cell_std=std(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    err_population_response(x,:,:)=(err_population_response(x,:,:)-cell_mean)./cell_std;
end
end
%if pcd_denoise==1
%component_touse=3;
%[coeff,score,latent] = pca(population_response);
%end
t=templateSVM('KernelFunction','linear');
correct_decode_label=[reshape(repmat([1:8],6,1),1,48),reshape(repmat([1:8],6,1),1,48)];

%generate svm model using correct data
for i=1:length(edge1)
    OBJ=fitcecoc(correct_population_response(:,:,i)',correct_decode_label,'Coding','onevsall','Learners',t);
    cvecoc = crossval(OBJ,'KFold',10);
  %  cvecoc = crossval(OBJ,'leaveout','on');
    Yhat = kfoldPredict(cvecoc);
    confuse_mat=confusionmat(correct_decode_label,Yhat);
    num_correct=diag(confuse_mat);
    correct_pertrue_class=num_correct./sum(confuse_mat,2);
    mean_correct_linearSVM=mean(correct_pertrue_class);
    time_perform_all(i)=mean_correct_linearSVM;
    full_model{i}=cvecoc;
end
temp1=1+(condition_index-1)*6;
temp2=[temp1;temp1+1;temp1+2;temp1+3;temp1+4;temp1+5];
temp3=reshape(temp2,1,6*length(condition_index));

for i=1:length(edge1)
    model_touse=full_model{i};
    temp_partition=model_touse.Partition; %find test trails for each fold
    test_include_trial=zeros(10,length(temp3));
    numcorrect=[];
    for p=1:10
    fold_test_index=find(training(temp_partition,p)==0);
    [C,ia,ib]=intersect(fold_test_index,temp3);
    test_include_trial(p,ib)=1;
    end

    for j=1:length(temp3)
    cv_touse=find(test_include_trial(:,j)==1);
    Yhat = predict(model_touse.Trained{cv_touse},correct_population_response(:,temp3(j),i)');
    numcorrect(j)=length(find((sac_correct_label(j)-Yhat')==0));
    end
 %   time_perform_correct(i)=sum(numcorrect)/length(temp3);
    total_trials=length(condition_index)*num_trial_percondition;
    time_perform_correct(i)=sum(numcorrect(randperm(length(temp3),total_trials)))/total_trials;
end
err_lookup=[5,6,7,8,1,2,3,4];

for i=1:length(edge1)
    model_touse=full_model{i};
    numcorrect=[];
    for j=1:10
    Yhat = predict(model_touse.Trained{j},err_population_response(:,:,i)');
    numcorrect(j)=length(find((err_lookup(sac_err_label)-Yhat')==0));
    end
    time_perform_err(i)=mean(numcorrect)/length(Yhat);
end
x_time=linspace(-1,6,length(edge1))+0.5*win_size;
end
