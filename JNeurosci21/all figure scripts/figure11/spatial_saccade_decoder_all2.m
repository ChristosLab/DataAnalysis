function [x_time,time_perform_correct]=spatial_saccade_decoder_all2(win_size,step_size,correctcell_index)
%this decoder decod both correct and err trials
%load('pass_selectivity_sac.mat');
%load('pass_strength_sac.mat');
%load('all_saccells.mat');
%win_size=0.4;step_size=0.1;
%correctcell_index=intersect(all_saccells,pass_selectivity_sac);
%correctcell_index=intersect(all_saccells,pass_strength_sac);

load('all_spatial_data.mat');
load('spatial_sacdirect.mat');
%load('mid_dorsal_informative_pre.mat');
informative_select=correctcell_index;
%informative_select([504,505,506])=[]; %use when include all cells
%informative_select([823,824])=[];
normalize=1;  %if normalize the data

%step_size=0.1;
%win_size=0.2;
edge1=[-1:step_size:6-win_size];   %6 second of time 
edge2=[-1+win_size:step_size:6];
correct_cell_id=[];
%correct_cueloc=[];
correct_ifmatch=[];
correct_bincount=[];
correct_sacdirect=[];

%store_num_trials=[];
correct_trial_count=0;

%data extraction
for i=1:length(informative_select) %loop cell
    cell_saccade=[all_trialdirect{informative_select(i)}]; %all saccade direction of the cell
    correct_cell_data=all_spatial_data(informative_select(i),:);  
    for j=1:8  %loop cuecalss
        temp_correctdata=correct_cell_data{j};%get physioogy data
        correct_trialnum=[temp_correctdata.trialnum]; %correct trial index
        correct_class_sacdirect=cell_saccade(correct_trialnum); %correct trial sac direction
        
        correct_cell_id=[correct_cell_id,i*ones(1,length(temp_correctdata))]; %store task information
    %    correct_cueloc=[correct_cueloc,j*ones(1,length(temp_correctdata))];
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
    end
end

%build pseudopopulation
for i=1:length(informative_select) %loop through cell
   % disp(i);
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
%end
if normalize==1
for x=1:size(correct_population_response)
    temp_data=squeeze(correct_population_response(x,:,:));
    cell_mean=mean(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    cell_std=std(reshape(temp_data,1,size(temp_data,2)*size(temp_data,1)));
    correct_population_response(x,:,:)=(correct_population_response(x,:,:)-cell_mean)./cell_std;
end
end
%if pcd_denoise==1
%component_touse=3;
%[coeff,score,latent] = pca(population_response);
%end
t=templateSVM('KernelFunction','linear');
%svm for each time point
decode_label=[reshape(repmat([1:8],6,1),1,48),reshape(repmat([1:8],6,1),1,48)];
for i=1:length(edge1)
   % disp(i);

    OBJ=fitcecoc(correct_population_response(:,:,i)',decode_label,'Coding','onevsall','Learners',t);
    cvecoc = crossval(OBJ,'KFold',10);
   % cvecoc = crossval(OBJ,'leaveout','on');
    Yhat = kfoldPredict(cvecoc);
    confuse_mat=confusionmat(decode_label,Yhat);
    num_correct=diag(confuse_mat);
    correct_pertrue_class=num_correct./sum(confuse_mat,2);
    mean_correct_linearSVM=mean(correct_pertrue_class);
    time_perform_correct(i)=mean_correct_linearSVM;
end

x_time=linspace(-1,6,length(edge1))+0.5*win_size;
end
