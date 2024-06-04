clear all; close all; clc;
%this version use all trials from both correct and error conditions to
%calculate base and for final plot

%%% reading all the mat files from the directory
BrainArea='PFC';
[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Rana Scripts\ODRdistVar Rana\ODRdistVar',BrainArea);


%%

neuron=0;
indices=[];
for i=1:length(raw)
    %load the matfiles
    fn_corr = [raw{i,1}, '_', num2str(raw{i,2}), '.mat'];
    cell_data = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);

    fn_err = [raw{i,1}, '_', num2str(raw{i,2}), '_err.mat'];
    cell_err_data =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);

    %calculate best class and sort classes
    max_class_corr = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    if max_class_corr(1) == 1
        Classes = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
    elseif max_class_corr(1) ==6
        Classes = [6 7 8 9 10 1 2 3 4 5 16 17 18 19 20 11 12 13 14 15];
    else
%         good_data(i)=0;
        continue;
    end
    corr_classes=0; err_classes=0;
    for j=Classes(1:5)
        try
            corr_classes=corr_classes+length(cell_data.MatData.class(j).ntr);
        catch
        end
        try
            err_classes=err_classes+length(cell_err_data.MatData.class(j).ntr);
        catch
        end
    end

    if corr_classes*err_classes>0
        corr_classes=0; err_classes=0;
        for j=Classes(6:10)
            try
                corr_classes=corr_classes+length(cell_data.MatData.class(j).ntr);
            catch
            end
            try
                err_classes=err_classes+length(cell_err_data.MatData.class(j).ntr);
            catch
            end
        end
        if corr_classes*err_classes>0
            corr_classes=0; err_classes=0;
            for j=Classes(11:15)
                try
                    corr_classes=corr_classes+length(cell_data.MatData.class(j).ntr);
                catch
                end
                try
                    err_classes=err_classes+length(cell_err_data.MatData.class(j).ntr);
                catch
                end
            end
            if corr_classes*err_classes>0
                corr_classes=0; err_classes=0;
                for j=Classes(16:20)
                    try
                        corr_classes=corr_classes+length(cell_data.MatData.class(j).ntr);
                    catch
                    end
                    try
                        err_classes=err_classes+length(cell_err_data.MatData.class(j).ntr);
                    catch
                    end
                end
                if corr_classes*err_classes>0
                    neuron=neuron+1;
                    indices=[indices i];
                end
            end
        end
    end
    clear cell_data cell_err_data;
end

clearvars -except raw indices neuron BrainArea

%%

base_time       = [3];   %time point to culculate pca base
good_data       = ones(1,length(indices));
winsize         = 0.5; %parameters for calcualting space
winstep         = 0.5;
num_timepoints  = 12;
winsize2        = 0.25; % parameters for trajectory
winstep2        = 0.1;
num_timepoints2 = 60;

firingrates     = nan(length(indices),20,num_timepoints, 20);%%%%%%%%%%%%%%%%%for caculating space
firingrates2    = nan(length(indices),20,num_timepoints2,20);%%%%%%%%%%%%%%%%%for caculating trajactory
firingrates_err = nan(length(indices),20,num_timepoints, 20);
firingrates2_err= nan(length(indices),20,num_timepoints2,20);

edge1=[-1:winstep:-1+(num_timepoints-1)*winstep];
edge2=[-1+winsize:winstep:-1+winsize+(num_timepoints-1)*winstep];
edge3=[-1:winstep2:-1+(num_timepoints2-1)*winstep2];
edge4=[-1+winsize2:winstep2:-1+winsize2+(num_timepoints2-1)*winstep2];%%%%%%%%%%%%calculate bin edge for space and trajactory

count_nofiring      =0;
count_extremefiring =0;
count_fewtrial      =0;
min_trialnum        =6;
analyse_trialnum    =6;
%%
for i=1:length(indices)

    %load the matfiles
    fn_corr = [raw{indices(i),1}, '_', num2str(raw{indices(i),2}), '.mat'];
    cell_data = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);

    fn_err = [raw{indices(i),1}, '_', num2str(raw{indices(i),2}), '_err.mat'];
    cell_err_data =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);

    %calculate best class and sort classes
    max_class_corr = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    if max_class_corr(1) == 1
        Classes = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
    elseif max_class_corr(1) ==6
        Classes = [6 7 8 9 10 1 2 3 4 5 16 17 18 19 20 11 12 13 14 15];
    else
        good_data(i)=0;
        continue;
    end

    if ~isempty(cell_data.MatData) && ~isempty(cell_err_data.MatData)
        for j=Classes
           
            %%%%%%%%%%%%%%% for correct %%%%%%%%%%%%%%%%
            try
                cue_class_trials(j)=length(cell_data.MatData.class(j).ntr);
                trialNum(i,j)=cue_class_trials(j);

                cueclass_cueon   =[cell_data.MatData.class(j).ntr.Cue_onT];
                cueclass_TS      ={cell_data.MatData.class(j).ntr.TS};

                for p=1:cue_class_trials(j)
                    temp_spiketime              =cueclass_TS{p}-cueclass_cueon(p);

                    for q=1:num_timepoints
                        firingrates(i,j,q,p)=length(find(temp_spiketime>edge1(q) & temp_spiketime<edge2(q)))*(1/winsize);
                    end
                    for q=1:num_timepoints2
                        %%%%%%%%%%%%%four dimensions are: cells, class,time points, trial
                        firingrates2(i,j,q,p)=length(find(temp_spiketime>edge3(q) & temp_spiketime<edge4(q)))*(1/winsize2);
                    end
                end
                clear temp_spiketime
            catch
            end

            %%%%%%%%%%%%%%%%%%%error trials%%%%%%%%%%%%%%
            try
                cue_class_errtrials(j)=length(cell_err_data.MatData.class(j).ntr);
                errtrialNum(i,j)=cue_class_errtrials(j);

                cueclass_err_cueon    =[cell_err_data.MatData.class(j).ntr.Cue_onT];
                cueclass_err_TS       ={cell_err_data.MatData.class(j).ntr.TS};

                for p=1:cue_class_errtrials(j)
                    temp_spiketime2   = cueclass_err_TS{p}-cueclass_err_cueon(p);
                    for q=1:num_timepoints
                        firingrates_err(i,j,q,p)=length(find(temp_spiketime2>edge1(q) & temp_spiketime2<edge2(q)))*(1/winsize);
                    end
                    for q=1:num_timepoints2
                        %%%%%%%%%%%%%five dimensions are: cells, cue location,decision status,time points, trial
                        firingrates2_err(i,j,q,p)=length(find(temp_spiketime2>edge3(q) & temp_spiketime2<edge4(q)))*(1/winsize2);
                    end
                end
                clear temp_spiketime2;

            catch
            end
        end
    else
        good_data(i)=0;
    end

    store_meanrates(i)=mean(firingrates(i,:,:,:),[2,3,4],'omitnan');

    if max(firingrates(i,:,:,:),[],'all')==0 || isnan(max(firingrates(i,:,:,:,:),[],'all'))  %%%%%%%%%%%%%filter bad cells
        good_data(i)=0;
        count_nofiring=count_nofiring+1;
    end

%     if store_meanrates(i)<1  %%%%%%%%%%%%%%%%%%%%%%filter extreme firingrates
%         good_data(i)=0;
%         count_extremefiring=count_extremefiring+1;
%     end
% 
%     try
%         if min(cue_class_trials)<min_trialnum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%filter trial number
%             good_data(i)=0;
%             count_fewtrial=count_fewtrial+1;
%         end
%     catch
%         good_data(i)=0;
%         count_fewtrial=count_fewtrial+1;
%     end
    clear cell_data cell_err_data;
end
        
%%
% %%%%%%%%%%%procedure for selecting random 6 trials everytime%%%%%%%%%%%%%%%%%%%%%%%%%%
% for x=1:size(firingrates,1)
%     for z=1:20
%         trial_num=trialNum(x,z);
%         if trial_num>=analyse_trialnum
%             temp_randselect=randperm(trial_num,analyse_trialnum);
%             firingrates2(x,z,:,1:analyse_trialnum)=firingrates2(x,z,:,temp_randselect);
%             firingrates(x,z,:,1:analyse_trialnum) =firingrates(x,z,:,temp_randselect);
%         end
%     end
% end
%%
bin_firingrates_R1=squeeze(firingrates2(find(good_data==1),1:10,:,:));%%%%%%%%%%%%%%%%for trajactory, only use match trials and selected number of trials
trial_fix=repmat(mean(squeeze(firingrates2(find(good_data==1),1:10,1:7,:)),3,'omitnan'),1,1,num_timepoints2,1);%calculate trial fix rate to align fix period
bin_firingrates_R1=bin_firingrates_R1-trial_fix;

bin_firingrates_R2=squeeze(firingrates2(find(good_data==1),11:20,:,:));%%%%%%%%%%%%%%%%for trajactory, only use match trials and selected number of trials
trial_fix=repmat(mean(squeeze(firingrates2(find(good_data==1),11:20,1:7,:)),3,'omitnan'),1,1,num_timepoints2,1);%calculate trial fix rate to align fix period
bin_firingrates_R2=bin_firingrates_R2-trial_fix;

bin_firingrates=cat(2,bin_firingrates_R1,bin_firingrates_R2);

%%
bin_firingrates_R1_best =nanmean(bin_firingrates(:,1:5,:,:),2);
bin_firingrates_R1_dia  =nanmean(bin_firingrates(:,6:10,:,:),2);
bin_firingrates_R2_best =nanmean(bin_firingrates(:,11:15,:,:),2);
bin_firingrates_R2_dia  =nanmean(bin_firingrates(:,16:20,:,:),2);

bin_firingrates=cat(2,bin_firingrates_R1_best,bin_firingrates_R1_dia,bin_firingrates_R2_best,bin_firingrates_R2_dia);

%%
bin_firingrates_R1_err=squeeze(firingrates2_err(find(good_data==1),1:10,:,:));%%%%%%%%%%%%%%%%for trajactory, only use match trials and selected number of trials
trial_fix=repmat(mean(squeeze(firingrates2_err(find(good_data==1),1:10,1:7,:)),3,'omitnan'),1,1,num_timepoints2,1);%calculate trial fix rate to align fix period
bin_firingrates_R1_err=bin_firingrates_R1_err-trial_fix;

bin_firingrates_R2_err=squeeze(firingrates2_err(find(good_data==1),11:20,:,:));%%%%%%%%%%%%%%%%for trajactory, only use match trials and selected number of trials
trial_fix=repmat(mean(squeeze(firingrates2_err(find(good_data==1),11:20,1:7,:)),3,'omitnan'),1,1,num_timepoints2,1);%calculate trial fix rate to align fix period
bin_firingrates_R2_err=bin_firingrates_R2_err-trial_fix;

bin_firingrates_err=cat(2,bin_firingrates_R1_err,bin_firingrates_R2_err);

%%
bin_firingrates_R1_best_err =nanmean(bin_firingrates_err(:,1:5,:,:),2);
bin_firingrates_R1_dia_err  =nanmean(bin_firingrates_err(:,6:10,:,:),2);
bin_firingrates_R2_best_err =nanmean(bin_firingrates_err(:,11:15,:,:),2);
bin_firingrates_R2_dia_err  =nanmean(bin_firingrates_err(:,16:20,:,:),2);

bin_firingrates_err=cat(2,bin_firingrates_R1_best_err,bin_firingrates_R1_dia_err,bin_firingrates_R2_best_err,bin_firingrates_R2_dia_err);

%% fix substracted
all_firingrates_R1=squeeze(nanmean(firingrates(find(good_data==1), 1:10,base_time,:),3))-squeeze(mean(firingrates(find(good_data==1), 1:10,1:2,:),3,'omitnan'));%%%%%%%%%%%%%%%for space calculate,use match trial only and selected number of trials
all_firingrates_R2=squeeze(nanmean(firingrates(find(good_data==1), 11:20,base_time,:),3))-squeeze(mean(firingrates(find(good_data==1),11:20,1:2,:),3,'omitnan'));%%%%%%%%%%%%%%%for space calculate,use match trial only and selected number of trials
all_firingrates=cat(2,all_firingrates_R1,all_firingrates_R2);

%%
all_firingrates_R1_best =nanmean(all_firingrates(:,1:5,:,:),2);
all_firingrates_R1_dia  =nanmean(all_firingrates(:,6:10,:,:),2);
all_firingrates_R2_best =nanmean(all_firingrates(:,11:15,:,:),2);
all_firingrates_R2_dia  =nanmean(all_firingrates(:,16:20,:,:),2);

all_firingrates=cat(2,all_firingrates_R1_best,all_firingrates_R1_dia,all_firingrates_R2_best,all_firingrates_R2_dia);

%%
% all_firingrates_err_R1=squeeze(nanmean(firingrates_err(find(good_data==1),1:10,base_time,:),3))-squeeze(mean(firingrates_err(find(good_data==1),1:10,1:2,:),3,'omitnan'));
% all_firingrates_err_R2=squeeze(nanmean(firingrates_err(find(good_data==1),1:10,base_time,:),3))-squeeze(mean(firingrates_err(find(good_data==1),1:10,1:2,:),3,'omitnan'));
% 
% all_firingrates_err=cat(2, all_firingrates_err_R1,all_firingrates_err_R2);

% %%
% all_firingrates_R1_best_err =nanmean(all_firingrates_err(:,1:5,:,:),2);
% all_firingrates_R1_dia_err  =nanmean(all_firingrates_err(:,6:10,:,:),2);
% all_firingrates_R2_best_err =nanmean(all_firingrates_err(:,11:15,:,:),2);
% all_firingrates_R2_dia_err  =nanmean(all_firingrates_err(:,15:20,:,:),2);
% 
% all_firingrates_err=cat(2,all_firingrates_R1_best_err,all_firingrates_R1_dia_err,all_firingrates_R2_best_err,all_firingrates_R2_dia_err);
%%
for randselect=1:100
    temp_indices = randperm(size(all_firingrates,1), round(size(all_firingrates,1)*0.75));

    all_firingrates2     = all_firingrates(temp_indices, :, :);
%     all_firingrates_err2 = all_firingrates_err_all(indices, :, :);
    bin_firingrates_err2 = bin_firingrates_err(temp_indices, :, :, :);
    bin_firingrates2     = bin_firingrates(temp_indices, :, :, :);

    numcell=size(all_firingrates2,1);
    for q=1:numcell
        test_select_err(q,:,:,:)=bin_firingrates_err2(q,:,:,:);
        train_resample(q,:,:)=all_firingrates2(q,:,:);%%%%%%%%%%%choose half trials to avoid bias to epch that used to calculate space
        test_all(q,:,:,:)=bin_firingrates2(q,:,:,:);
    end
    %%
    train_mean=mean(train_resample,3,'omitnan');
    %normalize
    train_PCAmat=normalize(train_mean(:,:)','zscore');  %%%%%%%%%%%%%%%%%%%train data was normalized
    train_PCAmat(isnan(train_PCAmat))=0;
    store_mean=mean(train_mean(:,:)');
    store_std=std(train_mean(:,:)');
    
    [coeff,score,latent,tsquared,explained,mu] = pca(train_PCAmat);
    
    ex_variance(1:size(explained,1))=explained;
    No_neurons=size(train_PCAmat,2);
    
    PC1=coeff(:,1);
    PC2=coeff(:,2);
    PC3=coeff(:,3);
    
    all_data = mean(test_all, 4,'omitnan');%%%%%%%%%%%%trial averaged data,cell x cueloc x time
    
    select_data_err = mean(test_select_err, 4,'omitnan');
    
    % normalize
    all_timepoints_mat=all_data-repmat(store_mean',1,size(all_data,2),size(all_data,3));
    all_timepoints_mat=all_timepoints_mat./repmat(store_std',1,size(all_data,2),size(all_data,3));
    all_timepoints_mat(isnan(all_timepoints_mat))=0;
    
    all_timepoints_mat_err=select_data_err-repmat(store_mean',1,size(select_data_err,2),size(select_data_err,3));
    all_timepoints_mat_err=all_timepoints_mat_err./repmat(store_std',1,size(select_data_err,2),size(select_data_err,3));
    all_timepoints_mat_err(isnan(all_timepoints_mat_err))=0;
    
    %%
    for i=1:4
        temp_data=squeeze(all_timepoints_mat(:,i,:));
        PC1_projection=PC1'*temp_data;
        PC2_projection=PC2'*temp_data;
        corr_PC1s(randselect,i,:)=PC1_projection;
        corr_PC2s(randselect,i,:)=PC2_projection;
    end
    
    %%
    for i=1:4
        temp_data=squeeze(all_timepoints_mat_err(:,i,:));
        PC1_projection=PC1'*temp_data;
        PC2_projection=PC2'*temp_data;
        err_PC1s(randselect,i,:)=PC1_projection;
        err_PC2s(randselect,i,:)=PC2_projection;
    end
    disp(No_neurons)
    clearvars -except No_neurons all_firingrates bin_firingrates bin_firingrates_err randselect corr_PC1s corr_PC2s err_PC2s err_PC1s raw indices neuron BrainArea ex_variance
end


%% distance all the way
for randselect=1:100
% 
%     corr_r1_best =[squeeze(corr_PC1s(randselect,1,:)) squeeze(corr_PC2s(randselect,1,:))];
%     corr_r1_dia  =[squeeze(corr_PC1s(randselect,2,:)) squeeze(corr_PC2s(randselect,2,:))];
%     corr_r2_best =[squeeze(corr_PC1s(randselect,3,:)) squeeze(corr_PC2s(randselect,3,:))];
%     corr_r2_dia  =[squeeze(corr_PC1s(randselect,4,:)) squeeze(corr_PC2s(randselect,4,:))];
    R1_best_corr_cd=[mean(corr_PC1s(randselect,1,16:30),3),mean(corr_PC2s(randselect,1,16:30),3)];
    R1_dia_corr_cd =[mean(corr_PC1s(randselect,2,16:30),3),mean(corr_PC2s(randselect,2,16:30),3)];
    R2_best_corr_cd=[mean(corr_PC1s(randselect,3,16:30),3),mean(corr_PC2s(randselect,3,16:30),3)];
    R2_dia_corr_cd =[mean(corr_PC1s(randselect,4,16:30),3),mean(corr_PC2s(randselect,4,16:30),3)];

    R1_best_err_cd=[mean(err_PC1s(randselect,1,16:30),3),mean(err_PC2s(randselect,1,16:30),3)];
    R1_dia_err_cd =[mean(err_PC1s(randselect,2,16:30),3),mean(err_PC2s(randselect,2,16:30),3)];
    R2_best_err_cd=[mean(err_PC1s(randselect,3,16:30),3),mean(err_PC2s(randselect,3,16:30),3)];
    R2_dia_err_cd =[mean(err_PC1s(randselect,4,16:30),3),mean(err_PC2s(randselect,4,16:30),3)];

    corr_r1_dist(randselect,:) = norm(R1_best_corr_cd-R1_dia_corr_cd);
    corr_r2_dist(randselect,:) = norm(R2_best_corr_cd-R2_dia_corr_cd);
    err_r1_dist(randselect,:)  = norm(R1_best_err_cd-R1_dia_err_cd);
    err_r2_dist(randselect,:)  = norm(R2_best_err_cd-R2_dia_err_cd);

% 
%     err_r1_best =[squeeze(err_PC1s(randselect,1,:)) squeeze(err_PC2s(randselect,1,:))];
%     err_r1_dia  =[squeeze(err_PC1s(randselect,2,:)) squeeze(err_PC2s(randselect,2,:))];
%     err_r2_best =[squeeze(err_PC1s(randselect,3,:)) squeeze(err_PC2s(randselect,3,:))];
%     err_r2_dia  =[squeeze(err_PC1s(randselect,4,:)) squeeze(err_PC2s(randselect,4,:))];
% 
%     Cue =[mean(corr_r1_best(11:15,:))]

%     corr_r1_dist(randselect,:) = sqrt(sum((corr_r1_best-corr_r1_dia).^2,2))';
%     corr_r2_dist(randselect,:) = sqrt(sum((corr_r2_best-corr_r2_dia).^2,2))';
%     err_r1_dist(randselect,:)  = sqrt(sum((err_r1_best-err_r1_dia).^2,2))';
%     err_r2_dist(randselect,:)  = sqrt(sum((err_r2_best-err_r2_dia).^2,2))';

   
end
%%
% Cue   =[mean(corr_r1_dist(:,11:15),2) mean(err_r1_dist(:,11:15),2) mean(corr_r2_dist(:,11:15),2) mean(err_r2_dist(:,11:15),2)];
% CD    =[mean(corr_r1_dist(:,16:30),2) mean(err_r1_dist(:,16:30),2) mean(corr_r2_dist(:,16:30),2) mean(err_r2_dist(:,16:30),2)];
% Sample=[mean(corr_r1_dist(:,31:35),2) mean(err_r1_dist(:,31:35),2) mean(corr_r2_dist(:,31:35),2) mean(err_r2_dist(:,31:35),2)];
% SD    =[mean(corr_r1_dist(:,36:50),2) mean(err_r1_dist(:,36:50),2) mean(corr_r2_dist(:,36:50),2) mean(err_r2_dist(:,36:50),2)];

CD    =[corr_r1_dist err_r1_dist corr_r2_dist err_r2_dist];

%% paired ttest
% [h,p_value,ci,stats] = ttest(Cue(:,1), Cue(:,2))
% [h,p_value,ci,stats] = ttest(Cue(:,3), Cue(:,4))
% 
[h,p_value,ci,stats] = ttest(CD(:,1), CD(:,2))
[h,p_value,ci,stats] = ttest(CD(:,3), CD(:,4))
% 
% [h,p_value,ci,stats] = ttest(Sample(:,1), Sample(:,2))
% [h,p_value,ci,stats] = ttest(Sample(:,3), Sample(:,4))
% 
% [h,p_value,ci,stats] = ttest(SD(:,1), SD(:,2))
% [h,p_value,ci,stats] = ttest(SD(:,3), SD(:,4))



%%


CD(:);

figure
plot(CD(:,1)); hold on; plot(CD(:,2));

figure
plot(CD(:,3)); hold on; plot(CD(:,4));





