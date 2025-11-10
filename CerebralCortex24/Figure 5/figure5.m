clear all 
close all

load('test_data.mat');  %analyse pre or post training
load('test_info.mat'); 
mode=1; %mode1 decode cue, mode2 decode sample, mode3 decode match/nonmatch 
first_point_delay = [];
% FOR D1=0

for r=1:3
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
%5 cross
accuracy_mat=mean(all_results_5cross);
x=linspace(-1.1,5,length(y));%57 
xq=linspace(-1.1,5,400);
interp_mean_d0_5cross=smooth(interp1(x,accuracy_mat,xq),10);
%10 cross
accuracy_mat=mean(all_results_10cross);
x=linspace(-1.1,5,length(y));%57 
xq=linspace(-1.1,5,400);
interp_mean_d0_10cross=smooth(interp1(x,accuracy_mat,xq),10);
figure;
subplot(2,2,1);
plot(xq, interp_mean_d0_5cross, 'r', 'LineWidth', 2);
xlim([-1,5]);
ylabel('Accuracy','fontweight','bold');
xlabel('Time(s)','fontweight','bold');
title('Non-persistent Neurons','fontweight','bold');

subplot(2,2,3);
plot(xq, interp_mean_d0_10cross, 'r', 'LineWidth', 2);
xlim([-1,5]);
ylabel('Accuracy','fontweight','bold');
xlable('Time(s)','fontweight','bold');

% FOR D1=1

for r=1:3
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

%5cross
accuracy_mat=mean(all_results_5cross);
x=linspace(-1.1,5,length(y));%57 
xq=linspace(-1.1,5,400);
interp_mean_d1_5cross=smooth(interp1(x,accuracy_mat,xq),10);
%10 cross
accuracy_mat=mean(all_results_10cross);
x=linspace(-1.1,5,length(y));%57 
xq=linspace(-1.1,5,400);
interp_mean_d1_10cross=smooth(interp1(x,accuracy_mat,xq),10);


save d1_5cross.mat all_results_5cross
save d1_10cross.mat all_results_10cross

figure;
subplot(2,2,1);
plot(xq, interp_mean_d1_5cross, 'b', 'LineWidth', 2);
xlim([-1,5]);
ylabel('Accuracy','fontweight','bold');
xlabel('Time(s)','fontweight','bold');
title('Persistent Neurons','fontweight','bold');

subplot(2,2,3);
plot(xq, interp_mean_d1_10cross, 'b', 'LineWidth', 2);
xlim([-1,5]);
ylabel('Accuracy','fontweight','bold');
xlable('Time(s)','fontweight','bold');