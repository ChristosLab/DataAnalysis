clear all; close all; clc;
%this version use all trials from both correct and error conditions to
%calculate base and for final plot

%%% reading all the mat files from the directory
BrainArea='PPC';
[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Rana Scripts\ODRdistVar Rana\ODRdistVar',BrainArea);

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

                cueclass_cueon = [cell_data.MatData.class(j).ntr.Cue_onT];
                cueclass_TS    = {cell_data.MatData.class(j).ntr.TS};

                for p=1:cue_class_trials(j)
                    temp_spiketime = cueclass_TS{p}-cueclass_cueon(p);

                    for q=1:num_timepoints
                        firingrates(i,j,q,p) = length(find(temp_spiketime>edge1(q) & temp_spiketime<edge2(q)))*(1/winsize);
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

    if max(firingrates(i,:,:),[],'all')==0 || isnan(max(firingrates(i,:,:,:),[],'all'))  %%%%%%%%%%%%%filter bad cells
        good_data(i)=0;
        count_nofiring=count_nofiring+1;
    end
% 
%     if store_meanrates(i)<1  %%%%%%%%%%%%%%%%%%%%%%filter extreme firingrates
%         good_data(i)=0;
%         count_extremefiring=count_extremefiring+1;
%     end
% 
    try
        if min(cue_class_trials)<min_trialnum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%filter trial number
            good_data(i)=0;
            count_fewtrial=count_fewtrial+1;
        end
    catch
        good_data(i)=0;
        count_fewtrial=count_fewtrial+1;
    end
    clear cell_data cell_err_data;
end
        
%% procedure for selecting random 6 trials everytime
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
all_firingrates_R1_best =nanmean(all_firingrates(:,1:5,:),2);
all_firingrates_R1_dia  =nanmean(all_firingrates(:,6:10,:),2);
all_firingrates_R2_best =nanmean(all_firingrates(:,11:15,:),2);
all_firingrates_R2_dia  =nanmean(all_firingrates(:,16:20,:),2);

all_firingrates=cat(2,all_firingrates_R1_best,all_firingrates_R1_dia,all_firingrates_R2_best,all_firingrates_R2_dia);

%%
% all_firingrates_err_R1=squeeze(nanmean(firingrates_err(find(good_data==1),1:10,base_time,:),3))-squeeze(mean(firingrates_err(find(good_data==1),1:10,1:2,:),3,'omitnan'));
% all_firingrates_err_R2=squeeze(nanmean(firingrates_err(find(good_data==1),1:10,base_time,:),3))-squeeze(mean(firingrates_err(find(good_data==1),1:10,1:2,:),3,'omitnan'));
% 
% all_firingrates_err=cat(2, all_firingrates_err_R1,all_firingrates_err_R2);
% 
% %%
% all_firingrates_R1_best_err =nanmean(all_firingrates_err(:,1:5,:,:),2);
% all_firingrates_R1_dia_err  =nanmean(all_firingrates_err(:,6:10,:,:),2);
% all_firingrates_R2_best_err =nanmean(all_firingrates_err(:,11:15,:,:),2);
% all_firingrates_R2_dia_err  =nanmean(all_firingrates_err(:,15:20,:,:),2);
% 
% all_firingrates_err=cat(2,all_firingrates_R1_best_err,all_firingrates_R1_dia_err,all_firingrates_R2_best_err,all_firingrates_R2_dia_err);
%%
numcell=size(all_firingrates,1);
for q=1:numcell
    test_select_err(q,:,:,:)=bin_firingrates_err(q,:,:,:);
    train_resample(q,:,:)=all_firingrates(q,:,:);%%%%%%%%%%%choose half trials to avoid bias to epoch that used to calculate space
    test_all(q,:,:,:)=bin_firingrates(q,:,:,:);
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
    corr_PC1s(i,:)=PC1_projection;
    corr_PC2s(i,:)=PC2_projection;
end

%%
for i=1:4
    temp_data=squeeze(all_timepoints_mat_err(:,i,:));
    PC1_projection=PC1'*temp_data;
    PC2_projection=PC2'*temp_data;
    err_PC1s(i,:)=PC1_projection;
    err_PC2s(i,:)=PC2_projection;
end
%%
color_mat=[...
    0, 0, 0.55;             0.25, 0.41, 0.88;...
    0.39, 0.58, 0.93;       0.27, 0.51, 0.71;...
    0.68, 0.85, 0.90;       1, 0.5, 0.5;...
    1, 0.25, 0.25;          1, 0.41, 0.38;...
    0.88, 0.07, 0.37;       0.86, 0.08, 0.24;...
    0, 0, 0.55;             0.25, 0.41, 0.88;...
    0.39, 0.58, 0.93;       0.27, 0.51, 0.71;...
    0.68, 0.85, 0.90;       1, 0.5, 0.5;...
    1, 0.25, 0.25;          1, 0.41, 0.38;...
    0.88, 0.07, 0.37;       0.86, 0.08, 0.24;...
    ];   %using 8 colors to represent 16 conditions

%% averaged correct and error trajectories
time_axis=linspace(-1,5,num_timepoints2);


figure
hlegend1=[];
hlegend2=[];

subplot(1,2,1);
p1=plot3(smooth(corr_PC1s(1,:),10),smooth(corr_PC2s(1,:),10),smooth(time_axis),'color',color_mat(1,:),'LineWidth',2);
hlegend1=[hlegend1 p1];

hold on
p2=plot3(smooth(corr_PC1s(2,:),10),smooth(corr_PC2s(2,:),10),smooth(time_axis),'color',color_mat(10,:),'LineWidth',2);
hlegend1=[hlegend1 p2];

hold on
p3=plot3(smooth(err_PC1s(1,:),10),smooth(err_PC2s(1,:),10),smooth(time_axis),':','color',color_mat(1,:),'LineWidth',2);
hlegend1=[hlegend1 p3];

hold on
p4=plot3(smooth(err_PC1s(2,:),10),smooth(err_PC2s(2,:),10),smooth(time_axis),':','color',color_mat(10,:),'LineWidth',2);
hlegend1=[hlegend1 p4];

xlabel('PC1')
ylabel('PC2')
zlabel('Time(s)')
title([BrainArea, ' Remember 1st trials'])
legend(hlegend1,'Best Correct','Diametric Correct','Best Error','Diametric Error')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hlegend1=[];
hlegend2=[];

subplot(1,2,2);
p1=plot3(smooth(corr_PC1s(3,:),10),smooth(corr_PC2s(3,:),10),smooth(time_axis),'color',color_mat(1,:),'LineWidth',2);
hlegend1=[hlegend1 p1];

hold on
p2=plot3(smooth(corr_PC1s(4,:),10),smooth(corr_PC2s(4,:),10),smooth(time_axis),'color',color_mat(10,:),'LineWidth',2);
hlegend1=[hlegend1 p2];

hold on
p3=plot3(smooth(err_PC1s(3,:),10),smooth(err_PC2s(3,:),10),smooth(time_axis),':','color',color_mat(1,:),'LineWidth',2);
hlegend1=[hlegend1 p3];

hold on
p4=plot3(smooth(err_PC1s(4,:),10),smooth(err_PC2s(4,:),10),smooth(time_axis),':','color',color_mat(10,:),'LineWidth',2);
hlegend1=[hlegend1 p4];

xlabel('PC1')
ylabel('PC2')
zlabel('Time(s)')
title([BrainArea, ' Remember 2nd trials'])
legend(hlegend1,'Best Correct','Diametric Correct','Best Error','Diametric Error')
hold off

if BrainArea=='PFC'
    sgtitle(sprintf('PFC Averaged Correct and Error trajectories, n=%d',No_neurons))
elseif BrainArea=='PPC'
    sgtitle(sprintf('PPC Averaged Correct and Error trajectories, n=%d',No_neurons))
end




%% Sample distance points Correct-Error

figure
% subplot 222

R1_best_corr_cd=[mean(corr_PC1s(1,16:30),2),mean(corr_PC2s(1,16:30),2)];
R1_dia_corr_cd =[mean(corr_PC1s(2,16:30),2),mean(corr_PC2s(2,16:30),2)];
R2_best_corr_cd=[mean(corr_PC1s(3,16:30),2),mean(corr_PC2s(3,16:30),2)];
R2_dia_corr_cd =[mean(corr_PC1s(4,16:30),2),mean(corr_PC2s(4,16:30),2)];

scatter(R1_best_corr_cd(1), R1_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1); hold on;
scatter(R1_dia_corr_cd(1), R1_dia_corr_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1); hold on;
scatter(R2_best_corr_cd(1), R2_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1); hold on;
scatter(R2_dia_corr_cd(1), R2_dia_corr_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1); hold on;

R1_best_err_cd=[mean(err_PC1s(1,16:30),2),mean(err_PC2s(1,16:30),2)];
R1_dia_err_cd =[mean(err_PC1s(2,16:30),2),mean(err_PC2s(2,16:30),2)];
R2_best_err_cd=[mean(err_PC1s(3,16:30),2),mean(err_PC2s(3,16:30),2)];
R2_dia_err_cd =[mean(err_PC1s(4,16:30),2),mean(err_PC2s(4,16:30),2)];

scatter(R1_best_err_cd(1), R1_best_err_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
scatter(R1_dia_err_cd(1), R1_dia_err_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
scatter(R2_best_err_cd(1), R2_best_err_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
scatter(R2_dia_err_cd(1), R2_dia_err_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;

line([R1_best_corr_cd(1), R1_dia_corr_cd(1)], [R1_best_corr_cd(2), R1_dia_corr_cd(2)], 'Color','k','LineStyle','-','LineWidth',0.5); hold on;
line([R1_best_err_cd(1),  R1_dia_err_cd(1)],  [R1_best_err_cd(2),  R1_dia_err_cd(2)],  'Color','k','LineStyle','--','LineWidth',0.5); hold on;
line([R2_best_corr_cd(1), R2_dia_corr_cd(1)], [R2_best_corr_cd(2), R2_dia_corr_cd(2)], 'Color','b','LineStyle','-','LineWidth',0.5); hold on;
line([R2_best_err_cd(1),  R2_dia_err_cd(1)],  [R2_best_err_cd(2),  R2_dia_err_cd(2)],  'Color','b','LineStyle','--','LineWidth',0.5); hold on;

xlabel('PC1')
ylabel('PC2')
% title([BrainArea,' cue delay'])
% title([BrainArea])
grid on
hold off

% measure angle
% v1 = [R1_best_corr_cd-R1_dia_corr_cd];
% v2 = [R2_best_corr_cd-R2_dia_corr_cd];
% 
% % Calculate the dot product and magnitudes
% dot_product = dot(v1, v2);
% magnitude_v1 = norm(v1);
% magnitude_v2 = norm(v2);
% 
% % Calculate the angle in radians and degrees
% angle_rad = acos(dot_product / (magnitude_v1 * magnitude_v2));
% angle_deg = rad2deg(angle_rad);
% 
% % Display the result
% disp(['The angle between the lines is: ' num2str(angle_deg) ' degrees']);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot 221
% 
% R1_best_corr_cd=[mean(corr_PC1s(1,11:15),2),mean(corr_PC2s(1,11:15),2)];
% R1_dia_corr_cd =[mean(corr_PC1s(2,11:15),2),mean(corr_PC2s(2,11:15),2)];
% R2_best_corr_cd=[mean(corr_PC1s(3,11:15),2),mean(corr_PC2s(3,11:15),2)];
% R2_dia_corr_cd =[mean(corr_PC1s(4,11:15),2),mean(corr_PC2s(4,11:15),2)];
% 
% scatter(R1_best_corr_cd(1), R1_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R1_dia_corr_cd(1), R1_dia_corr_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R2_best_corr_cd(1), R2_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1); hold on;
% scatter(R2_dia_corr_cd(1), R2_dia_corr_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1); hold on;
% 
% R1_best_err_cd=[mean(err_PC1s(1,11:15),2),mean(err_PC2s(1,11:15),2)];
% R1_dia_err_cd =[mean(err_PC1s(2,11:15),2),mean(err_PC2s(2,11:15),2)];
% R2_best_err_cd=[mean(err_PC1s(3,11:15),2),mean(err_PC2s(3,11:15),2)];
% R2_dia_err_cd =[mean(err_PC1s(4,11:15),2),mean(err_PC2s(4,11:15),2)];
% 
% scatter(R1_best_err_cd(1), R1_best_err_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% scatter(R1_dia_err_cd(1), R1_dia_err_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% scatter(R2_best_err_cd(1), R2_best_err_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% scatter(R2_dia_err_cd(1), R2_dia_err_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% 
% line([R1_best_corr_cd(1), R1_dia_corr_cd(1)], [R1_best_corr_cd(2), R1_dia_corr_cd(2)], 'Color','k','LineStyle','-','LineWidth',0.5); hold on;
% line([R1_best_err_cd(1),  R1_dia_err_cd(1)],  [R1_best_err_cd(2),  R1_dia_err_cd(2)],  'Color','k','LineStyle','--','LineWidth',0.5); hold on;
% line([R2_best_corr_cd(1), R2_dia_corr_cd(1)], [R2_best_corr_cd(2), R2_dia_corr_cd(2)], 'Color','b','LineStyle','-','LineWidth',0.5); hold on;
% line([R2_best_err_cd(1),  R2_dia_err_cd(1)],  [R2_best_err_cd(2),  R2_dia_err_cd(2)],  'Color','b','LineStyle','--','LineWidth',0.5); hold on;
% 
% xlabel('PC1')
% ylabel('PC2')
% title([BrainArea,' cue'])
% grid on
% hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot 223
% 
% R1_best_corr_cd=[mean(corr_PC1s(1,31:35),2),mean(corr_PC2s(1,31:35),2)];
% R1_dia_corr_cd =[mean(corr_PC1s(2,31:35),2),mean(corr_PC2s(2,31:35),2)];
% R2_best_corr_cd=[mean(corr_PC1s(3,31:35),2),mean(corr_PC2s(3,31:35),2)];
% R2_dia_corr_cd =[mean(corr_PC1s(4,31:35),2),mean(corr_PC2s(4,31:35),2)];
% 
% scatter(R1_best_corr_cd(1), R1_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R1_dia_corr_cd(1), R1_dia_corr_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R2_best_corr_cd(1), R2_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1); hold on;
% scatter(R2_dia_corr_cd(1), R2_dia_corr_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1); hold on;
% 
% R1_best_err_cd=[mean(err_PC1s(1,31:35),2),mean(err_PC2s(1,31:35),2)];
% R1_dia_err_cd =[mean(err_PC1s(2,31:35),2),mean(err_PC2s(2,31:35),2)];
% R2_best_err_cd=[mean(err_PC1s(3,31:35),2),mean(err_PC2s(3,31:35),2)];
% R2_dia_err_cd =[mean(err_PC1s(4,31:35),2),mean(err_PC2s(4,31:35),2)];
% 
% scatter(R1_best_err_cd(1), R1_best_err_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% scatter(R1_dia_err_cd(1), R1_dia_err_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% scatter(R2_best_err_cd(1), R2_best_err_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% scatter(R2_dia_err_cd(1), R2_dia_err_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% 
% line([R1_best_corr_cd(1), R1_dia_corr_cd(1)], [R1_best_corr_cd(2), R1_dia_corr_cd(2)], 'Color','k','LineStyle','-','LineWidth',0.5); hold on;
% line([R1_best_err_cd(1),  R1_dia_err_cd(1)],  [R1_best_err_cd(2),  R1_dia_err_cd(2)],  'Color','k','LineStyle','--','LineWidth',0.5); hold on;
% line([R2_best_corr_cd(1), R2_dia_corr_cd(1)], [R2_best_corr_cd(2), R2_dia_corr_cd(2)], 'Color','b','LineStyle','-','LineWidth',0.5); hold on;
% line([R2_best_err_cd(1),  R2_dia_err_cd(1)],  [R2_best_err_cd(2),  R2_dia_err_cd(2)],  'Color','b','LineStyle','--','LineWidth',0.5); hold on;
% 
% xlabel('PC1')
% ylabel('PC2')
% title([BrainArea,' sample'])
% grid on
% hold off
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot 224
% 
% R1_best_corr_cd=[mean(corr_PC1s(1,36:50),2),mean(corr_PC2s(1,36:50),2)];
% R1_dia_corr_cd =[mean(corr_PC1s(2,36:50),2),mean(corr_PC2s(2,36:50),2)];
% R2_best_corr_cd=[mean(corr_PC1s(3,36:50),2),mean(corr_PC2s(3,36:50),2)];
% R2_dia_corr_cd =[mean(corr_PC1s(4,36:50),2),mean(corr_PC2s(4,36:50),2)];
% 
% scatter(R1_best_corr_cd(1), R1_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R1_dia_corr_cd(1), R1_dia_corr_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R2_best_corr_cd(1), R2_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1); hold on;
% scatter(R2_dia_corr_cd(1), R2_dia_corr_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1); hold on;
% 
% R1_best_err_cd=[mean(err_PC1s(1,36:50),2),mean(err_PC2s(1,36:50),2)];
% R1_dia_err_cd =[mean(err_PC1s(2,36:50),2),mean(err_PC2s(2,36:50),2)];
% R2_best_err_cd=[mean(err_PC1s(3,36:50),2),mean(err_PC2s(3,36:50),2)];
% R2_dia_err_cd =[mean(err_PC1s(4,36:50),2),mean(err_PC2s(4,36:50),2)];
% 
% scatter(R1_best_err_cd(1), R1_best_err_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% scatter(R1_dia_err_cd(1), R1_dia_err_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% scatter(R2_best_err_cd(1), R2_best_err_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% scatter(R2_dia_err_cd(1), R2_dia_err_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1,'MarkerEdgeAlpha',.3); hold on;
% 
% line([R1_best_corr_cd(1), R1_dia_corr_cd(1)], [R1_best_corr_cd(2), R1_dia_corr_cd(2)], 'Color','k','LineStyle','-','LineWidth',0.5); hold on;
% line([R1_best_err_cd(1),  R1_dia_err_cd(1)],  [R1_best_err_cd(2),  R1_dia_err_cd(2)],  'Color','k','LineStyle','--','LineWidth',0.5); hold on;
% line([R2_best_corr_cd(1), R2_dia_corr_cd(1)], [R2_best_corr_cd(2), R2_dia_corr_cd(2)], 'Color','b','LineStyle','-','LineWidth',0.5); hold on;
% line([R2_best_err_cd(1),  R2_dia_err_cd(1)],  [R2_best_err_cd(2),  R2_dia_err_cd(2)],  'Color','b','LineStyle','--','LineWidth',0.5); hold on;
% 
% xlabel('PC1')
% ylabel('PC2')
% title([BrainArea,' sample delay'])
% grid on
% hold off


if BrainArea=='PFC'
    sgtitle(sprintf('PFC Distance Points, n=%d',No_neurons))
elseif BrainArea=='PPC'
    sgtitle(sprintf('PPC Distance Points, n=%d',No_neurons))
end



%%
% corr_r1_best =[squeeze(corr_PC1s(1,:))' squeeze(corr_PC2s(1,:))'];
% corr_r1_dia  =[squeeze(corr_PC1s(2,:))' squeeze(corr_PC2s(2,:))'];
% corr_r2_best =[squeeze(corr_PC1s(3,:))' squeeze(corr_PC2s(3,:))'];
% corr_r2_dia  =[squeeze(corr_PC1s(4,:))' squeeze(corr_PC2s(4,:))'];
% 
% err_r1_best =[squeeze(err_PC1s(1,:))' squeeze(err_PC2s(1,:))'];
% err_r1_dia  =[squeeze(err_PC1s(2,:))' squeeze(err_PC2s(2,:))'];
% err_r2_best =[squeeze(err_PC1s(3,:))' squeeze(err_PC2s(3,:))'];
% err_r2_dia  =[squeeze(err_PC1s(4,:))' squeeze(err_PC2s(4,:))'];
% 
% corr_r1_dist = sqrt(sum((corr_r1_best-corr_r1_dia).^2,2))';
% corr_r2_dist = sqrt(sum((corr_r2_best-corr_r2_dia).^2,2))';
% err_r1_dist = sqrt(sum((err_r1_best-err_r1_dia).^2,2))';
% err_r2_dist  = sqrt(sum((err_r2_best-err_r2_dia).^2,2))';
% 
% 
% %%
% 
% Cue   =[mean(corr_r1_dist(:,11:15),2) mean(err_r1_dist(:,11:15),2) mean(corr_r2_dist(:,11:15),2) mean(err_r2_dist(:,11:15),2)];
% CD    =[mean(corr_r1_dist(:,16:30),2) mean(err_r1_dist(:,16:30),2) mean(corr_r2_dist(:,16:30),2) mean(err_r2_dist(:,16:30),2)];
% Sample=[mean(corr_r1_dist(:,31:35),2) mean(err_r1_dist(:,31:35),2) mean(corr_r2_dist(:,31:35),2) mean(err_r2_dist(:,31:35),2)];
% SD    =[mean(corr_r1_dist(:,36:50),2) mean(err_r1_dist(:,36:50),2) mean(corr_r2_dist(:,36:50),2) mean(err_r2_dist(:,36:50),2)];
% 
% %%
% figure
% subplot 221
% bar(Cue)
% title([BrainArea, ' Cue'])
% subplot 222
% bar(CD)
% title([BrainArea, ' Cuedelay'])
% subplot 223
% bar(Sample)
% title([BrainArea, ' Sample'])
% subplot 224
% bar(SD)
% title([BrainArea, ' Sample delay'])
% 
% 
% 
% 
% 
% 
% 
%% plot distance all the way correct-error

% Corr_r1_best=[smooth(corr_PC1s(1,:),10), smooth(corr_PC2s(1,:),10)];
% Err_r1_best =[smooth(err_PC1s(1,:) ,10), smooth(err_PC2s(1,:) ,10)];
% 
% Corr_r1_dia =[smooth(corr_PC1s(2,:),10), smooth(corr_PC2s(2,:),10)];
% Err_r1_dia  =[smooth(err_PC1s(2,:) ,10), smooth(err_PC2s(2,:) ,10)];
% 
% figure
% subplot 121
% plot(smooth(time_axis), sqrt(sum((Corr_r1_best-Corr_r1_dia).^2 ,2)), 'k', 'LineWidth',2)
% hold on
% plot(smooth(time_axis), sqrt(sum((Err_r1_best-Err_r1_dia).^2 ,2)), 'k', 'LineWidth',2, 'LineStyle','--')
% 
% %adding color in between the outlined sections 
% patch1 = patch([0, 0, 0.5, 0.5], [0, 25, 25 ,0], 'r','EdgeColor','none'); 
% set(patch1, 'FaceAlpha', 0.15) 
% patch1 = patch([2, 2, 2.5, 2.5], [0, 25, 25 ,0], 'r','EdgeColor','none'); 
% set(patch1, 'FaceAlpha', 0.15) 
% patch1 = patch([4, 4, 5, 5], [0, 25, 25 ,0], 'b','EdgeColor','none'); 
% set(patch1, 'FaceAlpha', 0.15)
% set(gca, 'XLim',[-1 4.5])
% xlabel('Time(s)')
% ylabel('Trajectory distance')
% title([BrainArea, ' R1'])
% legend('Correct', 'Error')
% 
% 
% 
% Corr_r2_best=[smooth(corr_PC1s(3,:),10), smooth(corr_PC2s(3,:),10)];
% Err_r2_best =[smooth(err_PC1s(3,:) ,10), smooth(err_PC2s(3,:) ,10)];
% 
% Corr_r2_dia =[smooth(corr_PC1s(4,:),10), smooth(corr_PC2s(4,:),10)];
% Err_r2_dia  =[smooth(err_PC1s(4,:) ,10), smooth(err_PC2s(4,:) ,10)];
% 
% subplot 122
% plot(smooth(time_axis), sqrt(sum((Corr_r2_best-Corr_r2_dia).^2 ,2)), 'b', 'LineWidth',2)
% hold on
% plot(smooth(time_axis), sqrt(sum((Err_r2_best-Err_r2_dia).^2 ,2)), 'b', 'LineWidth',2, 'LineStyle','--')
% 
% %adding color in between the outlined sections 
% patch1 = patch([0, 0, 0.5, 0.5], [0, 25, 25 ,0], 'r','EdgeColor','none'); 
% set(patch1, 'FaceAlpha', 0.15) 
% patch1 = patch([2, 2, 2.5, 2.5], [0, 25, 25 ,0], 'r','EdgeColor','none'); 
% set(patch1, 'FaceAlpha', 0.15) 
% patch1 = patch([4, 4, 5, 5], [0, 25, 25 ,0], 'b','EdgeColor','none'); 
% set(patch1, 'FaceAlpha', 0.15)
% set(gca, 'XLim',[-1 4.5])
% xlabel('Time(s)')
% ylabel('Trajectory distance')
% title([BrainArea, ' R2'])
% legend('Correct', 'Error')
% sgtitle([BrainArea, ' Distance vs Time (Correct-Error)'])








