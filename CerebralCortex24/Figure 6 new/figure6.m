clear all 
close all

load('test_data.mat');  %analyse pre or post training
load('test_info.mat'); 
mode=1; %mode1 decode cue, mode2 decode sample, mode3 decode match/nonmatch 

%200 neurons, k = 10 for ALL decoders

%d1=0 graph
for r=1:3
    disp(r);
       
    %temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    %svm
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder(0.4,0.1,temp_index,mode);
    all_results_svm(r,:)=y;
    %DA
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=discriminant_decoder(0.4,0.1,temp_index,mode);
    all_results_DA(r,:)=y;
    %linear
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=linear_decoder(0.4,0.1,temp_index,mode);
    all_results_linear(r,:)=y;
    %kernel
   % temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
   % temp_index=temp_index(randperm(length(temp_index),200)); 
   % [x,y,crv_cell,population_response]=kernel_decoder(0.4,0.1,temp_index,mode);
   % all_results_kernel(r,:)=y;
    %classification trees
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_reponse]=trees_decoder(0.4,0.1,temp_index,mode);
    all_results_trees(r,:)=y;
    %naive bayes
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=naivebayes_decoder(0.4,0.1,temp_index,mode);
    all_results_NB(r,:)=y;
end

x=linspace(-1.1,5,length(y));%57 
xq=linspace(-1.1,5,400);

accuracy_mat=mean(all_results_svm);
interp_mean_SVM=smooth(interp1(x,accuracy_mat,xq),10);

accuracy_mat=mean(all_results_DA);
interp_mean_DA=smooth(interp1(x,accuracy_mat,xq),10);

accuracy_mat=mean(all_results_linear);
interp_mean_linear=smooth(interp1(x,accuracy_mat,xq),10);

%accuracy_mat=mean(all_results_kernel);
%interp_mean_kernel=smooth(interp1(x,accuracy_mat,xq),10);

accuracy_mat=mean(all_results_trees);
interp_mean_trees=smooth(interp1(x,accuracy_mat,xq),10);

accuracy_mat=mean(all_results_NB);
interp_mean_NB=smooth(interp1(x,accuracy_mat,xq),10);

all_d0_decoder=[interp_mean_SVM';interp_mean_DA';interp_mean_NB';interp_mean_trees'];

figure;
plot(xq, interp_mean_SVM);
hold on
plot(xq, interp_mean_DA );
%plot(xq, interp_mean_linear);
%plot(xq, interp_mean_kernel);
plot(xq, interp_mean_NB);
plot(xq, interp_mean_trees);
hold off

%d1=1 graph
for r=1:3
    disp(r);
      
    %svm
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder(0.4,0.1,temp_index,mode);
    all_results_svm(r,:)=y;
    %DA
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=discriminant_decoder(0.4,0.1,temp_index,mode);
    all_results_DA(r,:)=y;
    %linear
    %temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    %temp_index=temp_index(randperm(length(temp_index),200)); 
    %[x,y,crv_cell,population_response]=linear_decoder(0.4,0.1,temp_index,mode);
    %all_results_linear(r,:)=y;
    %kernel
    %temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    %temp_index=temp_index(randperm(length(temp_index),200)); 
    %[x,y,crv_cell,population_response]=kernel_decoder(0.4,0.1,temp_index,mode);
    %all_results_kernel(r,:)=y;
    %classification trees
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_reponse]=trees_decoder(0.4,0.1,temp_index,mode);
    all_results_trees(r,:)=y;
    %naive bayes
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=naivebayes_decoder(0.4,0.1,temp_index,mode);
    all_results_NB(r,:)=y;
end

x=linspace(-1.1,5,length(y));%57 
xq=linspace(-1.1,5,400);

accuracy_mat=mean(all_results_svm);
interp_mean_SVM=smooth(interp1(x,accuracy_mat,xq),10);

accuracy_mat=mean(all_results_DA);
interp_mean_DA=smooth(interp1(x,accuracy_mat,xq),10);

%accuracy_mat=mean(all_results_linear);
%interp_mean_linear=smooth(interp1(x,accuracy_mat,xq),10);

%accuracy_mat=mean(all_results_kernel);
%interp_mean_kernel=smooth(interp1(x,accuracy_mat,xq),10);

accuracy_mat=mean(all_results_trees);
interp_mean_trees=smooth(interp1(x,accuracy_mat,xq),10);

accuracy_mat=mean(all_results_NB);
interp_mean_NB=smooth(interp1(x,accuracy_mat,xq),10);

all_d1_decoder=[interp_mean_SVM';interp_mean_DA';interp_mean_NB';interp_mean_trees'];
figure;
plot(xq, interp_mean_SVM);
hold on
plot(xq, interp_mean_DA);
%plot(xq, interp_mean_linear);
%plot(xq, interp_mean_kernel);
plot(xq, interp_mean_NB);
plot(xq, interp_mean_trees);
hold off

save all_d0decoder_result.mat all_d0_decoder
save all_d1decoder_result.mat all_d1_decoder


%%%%%%%%%%%%%%%%%%%%%plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xq=linspace(-1.1,5,400);
load('all_d0decoder_result.mat');
interp_mean_SVM=all_d0_decoder(1,:);
interp_mean_DA=all_d0_decoder(2,:);
interp_mean_NB=all_d0_decoder(3,:);
interp_mean_trees=all_d0_decoder(4,:);
figure;
plot(xq, interp_mean_SVM,'r','LineWidth', 2);
hold on
plot(xq, interp_mean_DA,'--r' ,'LineWidth', 2);
%plot(xq, interp_mean_linear);
%plot(xq, interp_mean_kernel);
plot(xq, interp_mean_NB,':r','LineWidth', 2);
plot(xq, interp_mean_trees,'-.r','LineWidth', 2);
hold off
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
xlim([-1,5]);
legend('SVM','Discriminant Analysis','Native Bayes','Classification trees');
legend boxoff

load('all_d1decoder_result.mat');
interp_mean_SVM=all_d1_decoder(1,:);
interp_mean_DA=all_d1_decoder(2,:);
interp_mean_NB=all_d1_decoder(3,:);
interp_mean_trees=all_d1_decoder(4,:);
figure;
plot(xq, interp_mean_SVM,'b','LineWidth', 2);
hold on
plot(xq, interp_mean_DA,'--b','LineWidth', 2);
%plot(xq, interp_mean_linear);
%plot(xq, interp_mean_kernel);
plot(xq, interp_mean_NB,':b','LineWidth', 2);
plot(xq, interp_mean_trees,'-.b','LineWidth', 2);
hold off
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
xlim([-1,5]);
legend('SVM','Discriminant Analysis','Native Bayes','Classification trees');
legend boxoff