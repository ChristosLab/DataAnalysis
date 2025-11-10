clear all 
close all

load('test_data.mat');  %analyse pre or post training
load('test_info.mat'); 
mode=1; %mode1 decode cue, mode2 decode sample, mode3 decode match/nonmatch 
first_point_delay = [];
% FOR D1=0

for r=1:10
    disp(r);

    %200 neurons 
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_5cross(0.4,0.1,temp_index,mode);
    all_results_5cross(r,:)=y;
    
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_10cross(0.4,0.1,temp_index,mode);
    all_results_10cross(r,:)=y;

end
save d0_5cross.mat all_results_5cross
save d0_10cross.mat all_results_10cross
%load d0_5cross.mat
%load d0_10cross.mat
%5 cross
accuracy_mat=mean(all_results_5cross);
x=linspace(-1.1,5,length(accuracy_mat));%57 
xq=linspace(-1.1,5,400);
interp_mean_d0_5cross=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d0_5cross=smooth(interp1(x,std(all_results_5cross,0,1),xq),10);
%10 cross
accuracy_mat=mean(all_results_10cross);
x=linspace(-1.1,5,length(accuracy_mat));%57 
xq=linspace(-1.1,5,400);
interp_mean_d0_10cross=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d0_10cross=smooth(interp1(x,std(all_results_10cross,0,1),xq),10);
figure;
subplot(1,2,1);
hold on
r1=plot(xq, interp_mean_d0_5cross, 'r', 'LineWidth', 2);
curve1 = interp_mean_d0_5cross' + interp_std_d0_5cross';
curve2 = interp_mean_d0_5cross' - interp_std_d0_5cross';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
xlim([-1,5]);
ylabel('Accuracy','fontweight','bold');
xlabel('Time(s)','fontweight','bold');
title('5-fold cross validation','fontweight','bold');
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')

subplot(1,2,2);
hold on
r2=plot(xq, interp_mean_d0_10cross, 'r', 'LineWidth', 2);
curve1 = interp_mean_d0_10cross' + interp_std_d0_10cross';
curve2 = interp_mean_d0_10cross' - interp_std_d0_10cross';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
xlim([-1,5]);
ylabel('Accuracy','fontweight','bold');
xlabel('Time(s)','fontweight','bold');
title('10-fold cross validation','fontweight','bold');
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
% FOR D1=1

for r=1:10
    disp(r);
  
    %200 neurons 
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_5cross(0.4,0.1,temp_index,mode);
    all_results_5cross(r,:)=y;

    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_10cross(0.4,0.1,temp_index,mode);
    all_results_10cross(r,:)=y;
end
save d1_5cross.mat all_results_5cross
save d1_10cross.mat all_results_10cross
%load d1_5cross.mat
%load d1_10cross.mat

%5cross
accuracy_mat=mean(all_results_5cross);
x=linspace(-1.1,5,length(accuracy_mat));%57 
xq=linspace(-1.1,5,400);
interp_mean_d1_5cross=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d1_5cross=smooth(interp1(x,std(all_results_5cross,0,1),xq),10);
%10 cross
accuracy_mat=mean(all_results_10cross);
x=linspace(-1.1,5,length(accuracy_mat));%57 
xq=linspace(-1.1,5,400);
interp_mean_d1_10cross=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d1_10cross=smooth(interp1(x,std(all_results_10cross,0,1),xq),10);

subplot(1,2,1);
ylim([0,0.8]);
b1=plot(xq, interp_mean_d1_5cross, 'b', 'LineWidth', 2);
curve1 = interp_mean_d1_5cross' + interp_std_d1_5cross';
curve2 = interp_mean_d1_5cross' - interp_std_d1_5cross';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');
xlim([-1,5]);
ylabel('Accuracy','fontweight','bold');
xlabel('Time(s)','fontweight','bold');
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')

subplot(1,2,2);
ylim([0,0.8]);
b2=plot(xq, interp_mean_d1_10cross, 'b', 'LineWidth', 2);
curve1 = interp_mean_d1_10cross' + interp_std_d1_10cross';
curve2 = interp_mean_d1_10cross' - interp_std_d1_10cross';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');
xlim([-1,5]);
ylabel('Accuracy','fontweight','bold');
xlabel('Time(s)','fontweight','bold');
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
subplot(1,2,1)
legend([r1,b1],'Non-persistent','Persistent')
legend boxoff