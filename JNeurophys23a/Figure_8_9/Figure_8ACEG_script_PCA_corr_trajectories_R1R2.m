clear all; close all; clc;

%%% reading all the mat files from the directory
BrainArea='PPC';
[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Rana Scripts\ODRdistVar Rana\ODRdistVar',BrainArea);

%% for selecting neurons w at least 1 error trial in each of the four groups
% neuron=0;
% indices=[];
% for i=1:length(raw)
%     %load the matfiles
%     fn_corr = [raw{i,1}, '_', num2str(raw{i,2}), '.mat'];
%     cell_data = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
% 
%     fn_err = [raw{i,1}, '_', num2str(raw{i,2}), '_err.mat'];
%     cell_err_data =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);
% 
%     %calculate best class and sort classes
%     max_class_corr = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
%     if max_class_corr(1) == 1
%         Classes = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
%     elseif max_class_corr(1) ==6
%         Classes = [6 7 8 9 10 1 2 3 4 5 16 17 18 19 20 11 12 13 14 15];
%     else
% %         good_data(i)=0;
%         continue;
%     end
%     corr_classes=0; err_classes=0;
%     for j=Classes(1:5)
%         try
%             corr_classes=corr_classes+length(cell_data.MatData.class(j).ntr);
%         catch
%         end
%         try
%             err_classes=err_classes+length(cell_err_data.MatData.class(j).ntr);
%         catch
%         end
%     end
% 
%     if corr_classes*err_classes>0
%         corr_classes=0; err_classes=0;
%         for j=Classes(6:10)
%             try
%                 corr_classes=corr_classes+length(cell_data.MatData.class(j).ntr);
%             catch
%             end
%             try
%                 err_classes=err_classes+length(cell_err_data.MatData.class(j).ntr);
%             catch
%             end
%         end
%         if corr_classes*err_classes>0
%             corr_classes=0; err_classes=0;
%             for j=Classes(11:15)
%                 try
%                     corr_classes=corr_classes+length(cell_data.MatData.class(j).ntr);
%                 catch
%                 end
%                 try
%                     err_classes=err_classes+length(cell_err_data.MatData.class(j).ntr);
%                 catch
%                 end
%             end
%             if corr_classes*err_classes>0
%                 corr_classes=0; err_classes=0;
%                 for j=Classes(16:20)
%                     try
%                         corr_classes=corr_classes+length(cell_data.MatData.class(j).ntr);
%                     catch
%                     end
%                     try
%                         err_classes=err_classes+length(cell_err_data.MatData.class(j).ntr);
%                     catch
%                     end
%                 end
%                 if corr_classes*err_classes>0
%                     neuron=neuron+1;
%                     indices=[indices i];
%                 end
%             end
%         end
%     end
%     clear cell_data cell_err_data;
% end
% 
% clearvars -except raw indices neuron BrainArea

%%

base_time       = [3];   %time point to culculate pca base
good_data       = ones(1,length(raw));
winsize         = 0.5; %parameters for calcualting space
winstep         = 0.5;
num_timepoints  = 12;
winsize2        = 0.25; % parameters for trajectory
winstep2        = 0.1;
num_timepoints2 = 60;

firingrates     = nan(length(raw),20,num_timepoints, 20);%%%%%%%%%%%%%%%%%for caculating space
firingrates2    = nan(length(raw),20,num_timepoints2,20);%%%%%%%%%%%%%%%%%for caculating trajactory

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
for i=1:length(raw)

    %load the matfiles
    fn_corr = [raw{i,1}, '_', num2str(raw{i,2}), '.mat'];
    cell_data = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);

    %calculate best class and sort classes
    max_class_corr = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    if max_class_corr(1) == 1
        Classes=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
    elseif max_class_corr(1) ==6
        Classes = [6 7 8 9 10 1 2 3 4 5 16 17 18 19 20 11 12 13 14 15];
    else
        break;
    end

    if ~isempty(cell_data.MatData)
        for j=Classes
            %%%%%%%%%%%%%%% for correct %%%%%%%%%%%%%%%%
            try
                cue_class_trials(j) = length(cell_data.MatData.class(j).ntr);
                trialNum(i,j)       = cue_class_trials(j);

                cueclass_cueon = [cell_data.MatData.class(j).ntr.Cue_onT];
                cueclass_TS    = {cell_data.MatData.class(j).ntr.TS};

                for p_r1=1:cue_class_trials(j)
                    temp_spiketime = cueclass_TS{p_r1}-cueclass_cueon(p_r1);

                    for q=1:num_timepoints
                        firingrates(i,j,q,p_r1) = length(find(temp_spiketime>edge1(q) & temp_spiketime<edge2(q)))*(1/winsize);
                    end
                    for q=1:num_timepoints2
                        %%%%%%%%%%%%%four dimensions are: cells, class,time points, trial
                        firingrates2(i,j,q,p_r1)=length(find(temp_spiketime>edge3(q) & temp_spiketime<edge4(q)))*(1/winsize2);
                    end
                end
                clear temp_spiketime
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
%%
%%%%%%%%%%%procedure for selecting random 6 trials everytime%%%%%%%%%%%%%%%%%%%%%%%%%%
% for x=1:size(firingrates,1)
%     for z=1:20
%         trial_num=trialNum(x,z);
%         if trial_num>=analyse_trialnum
%             temp_randselect=randperm(trial_num,analyse_trialnum);
%             firingrates2(x,z,:,1:analyse_trialnum)=firingrates2(x,z,:,1:min_trialnum);
%             firingrates(x,z,:,1:analyse_trialnum) =firingrates(x,z,:,1:min_trialnum);
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
all_firingrates_R1=squeeze(nanmean(firingrates(find(good_data==1), 1:10,base_time,:),3))-squeeze(mean(firingrates(find(good_data==1), 1:10,1:2,:),3,'omitnan'));%%%%%%%%%%%%%%%for space calculate,use match trial only and selected number of trials
all_firingrates_R2=squeeze(nanmean(firingrates(find(good_data==1), 11:20,base_time,:),3))-squeeze(mean(firingrates(find(good_data==1),11:20,1:2,:),3,'omitnan'));%%%%%%%%%%%%%%%for space calculate,use match trial only and selected number of trials
all_firingrates=cat(2,all_firingrates_R1,all_firingrates_R2);


%%
numcell=size(all_firingrates,1);
for q=1:numcell
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

% normalize
all_timepoints_mat=all_data-repmat(store_mean',1,size(all_data,2),size(all_data,3));
all_timepoints_mat=all_timepoints_mat./repmat(store_std',1,size(all_data,2),size(all_data,3));
all_timepoints_mat(isnan(all_timepoints_mat))=0;
%%
for i=1:20
    temp_data=squeeze(all_timepoints_mat(:,i,:));
    PC1_projection=PC1'*temp_data;
    PC2_projection=PC2'*temp_data;
    corr_PC1s(i,:)=PC1_projection;
    corr_PC2s(i,:)=PC2_projection;
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
    ];   %using 10 colors to represent 20 conditions
% color_mat=color_mat;

%%
time_axis=linspace(-1,5,60);
p_r1=[];p_r2=[];
for i=1:20
    if i<=10
        subplot(1,2,1);
        p1=plot3(smooth(corr_PC1s(i,:),10),smooth(corr_PC2s(i,:),10),smooth(time_axis),'color',color_mat(i,:),'LineWidth',2);
        xlabel('PC1');ylabel('PC2');zlabel('Time(s)');title('Remember 1st');
        p_r1=[p_r1 p1];
    else
        subplot(1,2,2);
        p2=plot3(smooth(corr_PC1s(i,:),10),smooth(corr_PC2s(i,:),10),smooth(time_axis),'color',color_mat(i,:),'LineWidth',2);
        xlabel('PC1');ylabel('PC2');zlabel('Time(s)');title('Remember 2nd');
        p_r2=[p_r2 p2];
    end
    hold on
    
end

legend(p_r1,'Class 1','Class 2','Class 3','Class 4','Class 5','Class 6','Class 7','Class 8','Class 9','Class 10')
legend(p_r2,'Class 11','Class 12','Class 13','Class 14','Class 15','Class 16','Class 17','Class 18','Class 19','Class 20')
if BrainArea=='PFC'
    sgtitle(sprintf('PFC Correct trajectories, n=%d',No_neurons))
elseif BrainArea=='PPC'
    sgtitle(sprintf('PPC Correct trajectories, n=%d',No_neurons))
end
%%
% figure
% boundary_shrink=0;
% % Example data: Replace these with your actual data
% x = [smooth(corr_PC1s(1,23),10) smooth(corr_PC1s(2,23),10) smooth(corr_PC1s(3,23),10) smooth(corr_PC1s(4,23),10) smooth(corr_PC1s(5,23),10) ...
%     smooth(corr_PC1s(6,23),10) smooth(corr_PC1s(7,23),10) smooth(corr_PC1s(8,23),10) smooth(corr_PC1s(9,23),10) smooth(corr_PC1s(10,23),10)];  % X-coordinates of the data points
% y = [smooth(corr_PC2s(1,23),10) smooth(corr_PC2s(2,23),10) smooth(corr_PC2s(3,23),10) smooth(corr_PC2s(4,23),10) smooth(corr_PC2s(5,23),10) ...
%     smooth(corr_PC2s(6,23),10) smooth(corr_PC2s(7,23),10) smooth(corr_PC2s(8,23),10) smooth(corr_PC2s(9,23),10) smooth(corr_PC2s(10,23),10)];  % Y-coordinates of the data points
% z = [0 0 0 0 0 0 0 0 0 0];  % Z-coordinates of the data points
% 
% % Create a 3D scatter plot
% scatter3(x, y, z, 'filled', 'MarkerFaceColor', 'k'); % 'filled' to fill markers, 'b' for blue color
% 
% % Draw lines connecting the points
% hold on;
% % hold on;
% % temp_vec=[x' y'];
% % k=boundary(temp_vec(:,1),temp_vec(:,2),boundary_shrink);
% % plot(temp_vec(k,1),temp_vec(k,2),'-k','LineWidth',1);
% % hold on;
% 
% % Example data: Replace these with your actual data
% x = [smooth(corr_PC1s(11,23),10) smooth(corr_PC1s(12,23),10) smooth(corr_PC1s(13,23),10) smooth(corr_PC1s(14,23),10) smooth(corr_PC1s(15,23),10) ...
%     smooth(corr_PC1s(16,23),10) smooth(corr_PC1s(17,23),10) smooth(corr_PC1s(18,23),10) smooth(corr_PC1s(19,23),10) smooth(corr_PC1s(20,23),10)];  % X-coordinates of the data points
% y = [smooth(corr_PC2s(11,23),10) smooth(corr_PC2s(12,23),10) smooth(corr_PC2s(13,23),10) smooth(corr_PC2s(14,23),10) smooth(corr_PC2s(15,23),10) ...
%     smooth(corr_PC2s(16,23),10) smooth(corr_PC2s(17,23),10) smooth(corr_PC2s(18,23),10) smooth(corr_PC2s(19,23),10) smooth(corr_PC2s(20,23),10)];  % Y-coordinates of the data points
% % z = [0 0 0 0 0 0 0 0 0 0];  % Z-coordinates of the data points
% 
% % Create a 3D scatter plot
% % scatter3(x, y, z, 'filled', 'MarkerFaceColor', 'b'); % 'filled' to fill markers, 'b' for blue color
% scatter(x, y, 'filled', 'MarkerFaceColor', 'b'); % 'filled' to fill markers, 'b' for blue color
% 
% % Draw lines connecting the points
% hold on;
% % temp_vec=[x' y'];
% % k=boundary(temp_vec(:,1),temp_vec(:,2),boundary_shrink);
% % plot(temp_vec(k,1),temp_vec(k,2),'-b','LineWidth',1);
% 
% 
% % Optional: Customize the plot
% xlabel('PC1');
% ylabel('PC2');
% zlabel('Time(s)');
% title('3D Scatter Plot with Lines');
% 
% % Optional: Add grid and adjust view
% grid on;
% view(3); % Set the 3D view for better visualization



%%
% time_axis=linspace(-1,5,num_timepoints2);
% figure
% hlegend1=[];
% hlegend2=[];
% for i=1:20
%     temp_data=squeeze(all_timepoints_mat(:,i,:));
%     PC1_projection=PC1'*temp_data;
%     PC2_projection=PC2'*temp_data;
%     corr_PC1s(i,:)=PC1_projection;
%     corr_PC2s(i,:)=PC2_projection;
%     if i<=10
%         subplot(1,2,1);
%         p1=plot3(smooth(PC1_projection,10),smooth(PC2_projection,10),smooth(time_axis),'color',color_mat(i,:),'LineWidth',2);
%         hlegend1=[hlegend1 p1];
%         xlabel('PC1')
%         ylabel('PC2')
%         zlabel('Time(sec)')
%         title([BrainArea, ' Remember 1st Correct trials'])
%         hold on
%     else
%         subplot(1,2,2);
%         p2=plot3(smooth(PC1_projection,10),smooth(PC2_projection,10),smooth(time_axis),'color',color_mat(i,:),'LineWidth',2);
%         hlegend2=[hlegend2 p2];
%         xlabel('PC1')
%         ylabel('PC2')
%         zlabel('Time(sec)')
%         title([BrainArea, ' Remember 2nd Correct trials'])
%         hold on
%     end
% end
% legend(hlegend1,'Class1','Class2','Class3','Class4','Class5','Class6','Class7','Class8','Class9','Class10')
% legend(hlegend2,'Class11','Class12','Class13','Class14','Class15','Class16','Class17','Class18','Class19','Class20')
% 
% 
% 
% %% averaged correct and error trajectories
% 
% figure
% hlegend1=[];
% hlegend2=[];
% 
% subplot(1,2,1);
% p1=plot3(smooth(mean(corr_PC1s(1:5,:),1),10),smooth(mean(corr_PC2s(1:5,:),1),10),smooth(time_axis),'color',color_mat(1,:),'LineWidth',2);
% hlegend1=[hlegend1 p1];
% 
% hold on
% p2=plot3(smooth(mean(corr_PC1s(6:10,:),1),10),smooth(mean(corr_PC2s(6:10,:),1),10),smooth(time_axis),'color',color_mat(10,:),'LineWidth',2);
% hlegend1=[hlegend1 p2];
% 
% xlabel('PC1')
% ylabel('PC2')
% zlabel('Time(sec)')
% title([BrainArea, ' Remember 1st trials(Cue)'])
% legend(hlegend1,'Best Correct','Diametric Correct')
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hlegend1=[];
% hlegend2=[];
% 
% subplot(1,2,2);
% p1=plot3(smooth(mean(corr_PC1s(11:15,:),1),10),smooth(mean(corr_PC2s(11:15,:),1),10),smooth(time_axis),'color',color_mat(1,:),'LineWidth',2);
% hlegend1=[hlegend1 p1];
% 
% hold on
% p2=plot3(smooth(mean(corr_PC1s(16:20,:),1),10),smooth(mean(corr_PC2s(16:20,:),1),10),smooth(time_axis),'color',color_mat(10,:),'LineWidth',2);
% hlegend1=[hlegend1 p2];
% 
% xlabel('PC1')
% ylabel('PC2')
% zlabel('Time(sec)')
% title([BrainArea, ' Remember 2nd trials(Cue)'])
% legend(hlegend1,'Best Correct','Diametric Correct')
% hold off
% 
% sgtitle([BrainArea, ' Averaged Correct trajectories'])


% %% Sample distance points Correct-Error
% 
% figure
% subplot 222
% 
% Corr_points=[mean(corr_PC1s(:,15:30),2),mean(corr_PC2s(:,15:30),2)];
% 
% R1_best_corr_cd = mean(Corr_points(1:5, :),1);
% R1_dia_corr_cd = mean(Corr_points(6:10, :),1);
% R2_best_corr_cd = mean(Corr_points(11:15, :),1);
% R2_dia_corr_cd = mean(Corr_points(16:20, :),1);
% 
% scatter(R1_best_corr_cd(1), R1_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R1_dia_corr_cd(1), R1_dia_corr_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R2_best_corr_cd(1), R2_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1); hold on;
% scatter(R2_dia_corr_cd(1), R2_dia_corr_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1); hold on;
% 
% line([R1_best_corr_cd(1), R1_dia_corr_cd(1)], [R1_best_corr_cd(2), R1_dia_corr_cd(2)], 'Color','k','LineStyle','-','LineWidth',0.5); hold on;
% line([R2_best_corr_cd(1), R2_dia_corr_cd(1)], [R2_best_corr_cd(2), R2_dia_corr_cd(2)], 'Color','b','LineStyle','-','LineWidth',0.5); hold on;
% 
% xlabel('PC1')
% ylabel('PC2')
% title([BrainArea,' cue'])
% grid on
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot 221
% 
% Corr_points=[mean(corr_PC1s(:,11:15),2),mean(corr_PC2s(:,11:15),2)];
% 
% R1_best_corr_cd = mean(Corr_points(1:5, :),1);
% R1_dia_corr_cd = mean(Corr_points(6:10, :),1);
% R2_best_corr_cd = mean(Corr_points(11:15, :),1);
% R2_dia_corr_cd = mean(Corr_points(16:20, :),1);
% 
% scatter(R1_best_corr_cd(1), R1_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R1_dia_corr_cd(1), R1_dia_corr_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R2_best_corr_cd(1), R2_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1); hold on;
% scatter(R2_dia_corr_cd(1), R2_dia_corr_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1); hold on;
% 
% line([R1_best_corr_cd(1), R1_dia_corr_cd(1)], [R1_best_corr_cd(2), R1_dia_corr_cd(2)], 'Color','k','LineStyle','-','LineWidth',0.5); hold on;
% line([R2_best_corr_cd(1), R2_dia_corr_cd(1)], [R2_best_corr_cd(2), R2_dia_corr_cd(2)], 'Color','b','LineStyle','-','LineWidth',0.5); hold on;
% 
% xlabel('PC1')
% ylabel('PC2')
% title([BrainArea,' cue'])
% grid on
% hold off
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot 223
% 
% Corr_points=[mean(corr_PC1s(:,31:35),2),mean(corr_PC2s(:,31:35),2)];
% 
% R1_best_corr_cd = mean(Corr_points(1:5, :),1);
% R1_dia_corr_cd = mean(Corr_points(6:10, :),1);
% R2_best_corr_cd = mean(Corr_points(11:15, :),1);
% R2_dia_corr_cd = mean(Corr_points(16:20, :),1);
% 
% scatter(R1_best_corr_cd(1), R1_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R1_dia_corr_cd(1), R1_dia_corr_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R2_best_corr_cd(1), R2_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1); hold on;
% scatter(R2_dia_corr_cd(1), R2_dia_corr_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1); hold on;
% 
% line([R1_best_corr_cd(1), R1_dia_corr_cd(1)], [R1_best_corr_cd(2), R1_dia_corr_cd(2)], 'Color','k','LineStyle','-','LineWidth',0.5); hold on;
% line([R2_best_corr_cd(1), R2_dia_corr_cd(1)], [R2_best_corr_cd(2), R2_dia_corr_cd(2)], 'Color','b','LineStyle','-','LineWidth',0.5); hold on;
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
% Corr_points=[mean(corr_PC1s(:,36:50),2),mean(corr_PC2s(:,36:50),2)];
% 
% R1_best_corr_cd = mean(Corr_points(1:5, :),1);
% R1_dia_corr_cd = mean(Corr_points(6:10, :),1);
% R2_best_corr_cd = mean(Corr_points(11:15, :),1);
% R2_dia_corr_cd = mean(Corr_points(16:20, :),1);
% 
% scatter(R1_best_corr_cd(1), R1_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R1_dia_corr_cd(1), R1_dia_corr_cd(2),100,'d','MarkerEdgeColor','k','LineWidth',1); hold on;
% scatter(R2_best_corr_cd(1), R2_best_corr_cd(2),100,'pentagram','MarkerEdgeColor','b','LineWidth',1); hold on;
% scatter(R2_dia_corr_cd(1), R2_dia_corr_cd(2),100,'d','MarkerEdgeColor','b','LineWidth',1); hold on;
% 
% line([R1_best_corr_cd(1), R1_dia_corr_cd(1)], [R1_best_corr_cd(2), R1_dia_corr_cd(2)], 'Color','k','LineStyle','-','LineWidth',0.5); hold on;
% line([R2_best_corr_cd(1), R2_dia_corr_cd(1)], [R2_best_corr_cd(2), R2_dia_corr_cd(2)], 'Color','b','LineStyle','-','LineWidth',0.5); hold on;
% 
% xlabel('PC1')
% ylabel('PC2')
% title([BrainArea,' sample'])
% grid on
% hold off
% 
% sgtitle([BrainArea, ' Distance Points'])



%%
% subplot(1,2,1);
% x1=xlim;
% y1=ylim;
% zlim([.2,.23]);
% view(0,90);
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% subplot(1,2,2);
% x2=xlim;
% y2=ylim;
% subplot(1,2,1);
% xlim([min([x1(1),x2(1)]),max([x1(2),x2(2)])]);
% ylim([min([y1(1),y2(1)]),max([y1(2),y2(2)])]);
% subplot(1,2,2);
% xlim([min([x1(1),x2(1)]),max([x1(2),x2(2)])]);
% ylim([min([y1(1),y2(1)]),max([y1(2),y2(2)])]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%plot temporal change
% figure;
% trialplot_list=[1,5];
% trial_slide=[1,2;2,3;3,4;4,5;5,6];
% % trial_slide=[1,2;2,3];
% trialorder_color_mat=[255 0 0;255 49 150;255 120 120;255 180 100;255 255 0]/255;
% for i=1:size(trial_slide,1)
%     tempdata_slide=bin_firingrates(:,:,:,trial_slide(i,:));
%     tempaverage_slide=nanmean(tempdata_slide, 4);
%     if if_normalize==1
%     %temptarg_mat=reshape(normalize(tempaverage_slide(:,:),2,'zscore'),size(tempaverage_slide));
%     temptarg_mat=(tempaverage_slide(:,:)-repmat(store_mean',1,size(tempaverage_slide(:,:),2)))./repmat(store_std',1,size(tempaverage_slide(:,:),2));
%     temptarg_mat=reshape(temptarg_mat,size(tempaverage_slide));
%     temptarg_mat(isinf(temptarg_mat)|isnan(temptarg_mat)) = 0;
%     else
%     temptarg_mat=tempaverage_slide;
%     end
%     for j=1:length(trialplot_list)
%        temp_data=squeeze(temptarg_mat(:,trialplot_list(j),:));
%        PC1_trialprojection(i,j,:)=PC1'*temp_data;
%        PC2_trialprojection(i,j,:)=PC2'*temp_data;
%     end
% end
%    
% %mean_PC1_trialprojection=mean(squeeze(mean(PC1_trialprojection)));
% %mean_PC2_trialprojection=mean(squeeze(mean(PC2_trialprojection))); 
% for i=1:size(trial_slide,1)
%    for j=1:length(trialplot_list)
%        %align fix period
%       %temp_PC1=squeeze(PC1_trialprojection(i,j,:));
%       %temp_PC2=squeeze(PC2_trialprojection(i,j,:));
%       %base_PC1=mean(temp_PC1(1:7));
%       %base_PC2=mean(temp_PC2(1:7));
%       %plot3(smooth(temp_PC1-base_PC1),smooth(temp_PC2-base_PC2),smooth(time_axis(:)),'color',trialorder_color_mat(i,:),'LineWidth',1);
%        plot3(smooth(squeeze(PC1_trialprojection(i,j,:)),10),smooth(squeeze(PC2_trialprojection(i,j,:)),10),smooth(time_axis(:),10),'color',trialorder_color_mat(i,:),'LineWidth',1.5);
%        hold on;
%    end
% end















































