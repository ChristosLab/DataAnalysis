clear all 
close all

load('test_data.mat');  %analyse pre or post training
load('test_info.mat'); 
mode=1; %mode1 decode cue, mode2 decode sample, mode3 decode match/nonmatch 
first_point_delay = [];
% FOR D1=0

for r=1:10
    disp(r);
    %select cell with delay activity   
    %temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    
    %50 neurons
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),50)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder(0.4,0.1,temp_index,mode);
    all_results50_d0(r,:)=y;
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),50)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_shuffle(0.4,0.1,temp_index,mode);
    all_results50_d0_shuffle(r,:)=y;
    
    %100 neurons
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),100)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder(0.4,0.1,temp_index,mode);
    all_results100_d0(r,:)=y;
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),100)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_shuffle(0.4,0.1,temp_index,mode);
    all_results100_d0_shuffle(r,:)=y;
    
    %200 neurons 
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder(0.4,0.1,temp_index,mode);
    all_results200_d0(r,:)=y;
    temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_shuffle(0.4,0.1,temp_index,mode);
    all_results200_d0_shuffle(r,:)=y;

end
%load d0_50cell.mat 
%load d0_100cell.mat 
%load d0_200cell.mat 
%load d0_50cell_shuffle.mat 
%load d0_100cell_shuffle.mat 
%load d0_200cell_shuffle.mat 

%50
accuracy_mat=mean(all_results50_d0);
accuracy_mat_shuffle=mean(all_results50_d0_shuffle);
x=linspace(-1.1,5,62);%57 
xq=linspace(-1.1,5,400);
interp_mean_d0_50cell=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d0_50cell=smooth(interp1(x,std(all_results50_d0,0,1),xq),10);
interp_mean_d0_50cell_shuffle=smooth(interp1(x,accuracy_mat_shuffle,xq),10);
interp_std_d0_50cell_shuffle=smooth(interp1(x,std(all_results50_d0_shuffle,0,1),xq),10);
first_point_delay(1) = interp_mean_d0_50cell(132);
%100
accuracy_mat=mean(all_results100_d0);
accuracy_mat_shuffle=mean(all_results100_d0_shuffle);
x=linspace(-1.1,5,62);%57 
xq=linspace(-1.1,5,400);
interp_mean_d0_100cell=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d0_100cell=smooth(interp1(x,std(all_results100_d0,0,1),xq),10);
interp_mean_d0_100cell_shuffle=smooth(interp1(x,accuracy_mat_shuffle,xq),10);
interp_std_d0_100cell_shuffle=smooth(interp1(x,std(all_results100_d0_shuffle,0,1),xq),10);
first_point_delay(2) = interp_mean_d0_100cell(132);
%200
accuracy_mat=mean(all_results200_d0);
accuracy_mat_shuffle=mean(all_results200_d0_shuffle);
x=linspace(-1.1,5,62);%57 
xq=linspace(-1.1,5,400);
interp_mean_d0_200cell=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d0_200cell=smooth(interp1(x,std(all_results200_d0,0,1),xq),10);
interp_mean_d0_200cell_shuffle=smooth(interp1(x,accuracy_mat_shuffle,xq),10);
interp_std_d0_200cell_shuffle=smooth(interp1(x,std(all_results200_d0_shuffle,0,1),xq),10);
first_point_delay(3) = interp_mean_d0_200cell(132);

save d0_50cell.mat all_results50_d0
save d0_100cell.mat all_results100_d0
save d0_200cell.mat all_results200_d0
save d0_50cell_shuffle.mat all_results50_d0_shuffle
save d0_100cell_shuffle.mat all_results100_d0_shuffle
save d0_200cell_shuffle.mat all_results200_d0_shuffle

% FOR D1=1

for r=1:10
    disp(r);
    %select cell with delay activity   
    
    
    %50 neurons
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),50)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder(0.4,0.1,temp_index,mode);
    all_results50_d1(r,:)=y;
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),50)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_shuffle(0.4,0.1,temp_index,mode);
    all_results50_d1_shuffle(r,:)=y;
    
    %100 neurons
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),100)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder(0.4,0.1,temp_index,mode);
    all_results100_d1(r,:)=y;
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),100)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_shuffle(0.4,0.1,temp_index,mode);
    all_results100_d1_shuffle(r,:)=y;

    %200 neurons 
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder(0.4,0.1,temp_index,mode);
    all_results200_d1(r,:)=y;
    temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
    temp_index=temp_index(randperm(length(temp_index),200)); 
    [x,y,crv_cell,population_response]=normalize_SVM_decoder_shuffle(0.4,0.1,temp_index,mode);
    all_results200_d1_shuffle(r,:)=y;

end

%load d1_50cell.mat 
%load d1_100cell.mat 
%load d1_200cell.mat 
%load d1_50cell_shuffle.mat 
%load d1_100cell_shuffle.mat 
%load d1_200cell_shuffle.mat 

%50
accuracy_mat=mean(all_results50_d1);
accuracy_mat_shuffle=mean(all_results50_d1_shuffle);
x=linspace(-1.1,5,length(x));%57 
xq=linspace(-1.1,5,400);
interp_mean_d1_50cell=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d1_50cell=smooth(interp1(x,std(all_results50_d1,0,1),xq),10);
interp_mean_d1_50cell_shuffle=smooth(interp1(x,accuracy_mat_shuffle,xq),10);
interp_std_d1_50cell_shuffle=smooth(interp1(x,std(all_results50_d1_shuffle,0,1),xq),10);
first_point_delay(2,1) = interp_mean_d1_50cell(132);
%100
accuracy_mat=mean(all_results100_d1);
accuracy_mat_shuffle=mean(all_results100_d1_shuffle);
x=linspace(-1.1,5,length(x));%57 
xq=linspace(-1.1,5,400);
interp_mean_d1_100cell=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d1_100cell=smooth(interp1(x,std(all_results100_d1,0,1),xq),10);
interp_mean_d1_100cell_shuffle=smooth(interp1(x,accuracy_mat_shuffle,xq),10);
interp_std_d1_100cell_shuffle=smooth(interp1(x,std(all_results100_d1_shuffle,0,1),xq),10);
first_point_delay(2,2) = interp_mean_d1_100cell(132);
%200
accuracy_mat=mean(all_results200_d1);
accuracy_mat_shuffle=mean(all_results200_d1_shuffle);
x=linspace(-1.1,5,length(x));%57 
xq=linspace(-1.1,5,400);
interp_mean_d1_200cell=smooth(interp1(x,accuracy_mat,xq),10);
interp_std_d1_200cell=smooth(interp1(x,std(all_results200_d1,0,1),xq),10);
interp_mean_d1_200cell_shuffle=smooth(interp1(x,accuracy_mat_shuffle,xq),10);
interp_std_d1_200cell_shuffle=smooth(interp1(x,std(all_results200_d1_shuffle,0,1),xq),10);
first_point_delay(2,3) = interp_mean_d1_200cell(132);

save d1_50cell.mat all_results50_d1
save d1_100cell.mat all_results100_d1
save d1_200cell.mat all_results200_d1
save d1_50cell_shuffle.mat all_results50_d1_shuffle
save d1_100cell_shuffle.mat all_results100_d1_shuffle
save d1_200cell_shuffle.mat all_results200_d1_shuffle

%figure4A
figure;
hold on;
r1=plot(xq,interp_mean_d0_50cell,'r','LineWidth',3);
b1=plot(xq,interp_mean_d1_50cell,'b','LineWidth',3);
curve1 = interp_mean_d0_50cell' + interp_std_d0_50cell';
curve2 = interp_mean_d0_50cell' - interp_std_d0_50cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
curve1 = interp_mean_d1_50cell' + interp_std_d1_50cell';
curve2 = interp_mean_d1_50cell' - interp_std_d1_50cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');
xlim([-1,5]);
yticks([0:0.2:0.8]);
ylim([0,0.8]);
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
legend([r1,b1],{'Non-persistent','Persistent'});
legend('boxoff');
title('Decoding with 50 Cells');

%figure4b
figure;
hold on;
r2=plot(xq,interp_mean_d0_100cell,'r','LineWidth',3);
b2=plot(xq,interp_mean_d1_100cell,'b','LineWidth',3);
curve1 = interp_mean_d0_100cell' + interp_std_d0_100cell';
curve2 = interp_mean_d0_100cell' - interp_std_d0_100cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
curve1 = interp_mean_d1_100cell' + interp_std_d1_100cell';
curve2 = interp_mean_d1_100cell' - interp_std_d1_100cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');
xlim([-1,5]);
yticks([0:0.2:0.8]);
ylim([0,0.8]);
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
legend([r2,b2],{'Non-persistent','Persistent'});
legend('boxoff');
title('Decoding with 100 Cells');

%figure4c
figure;
hold on;
r3=plot(xq,interp_mean_d0_200cell,'r','LineWidth',3);
b3=plot(xq,interp_mean_d1_200cell,'b','LineWidth',3);
curve1 = interp_mean_d0_200cell' + interp_std_d0_200cell';
curve2 = interp_mean_d0_200cell' - interp_std_d0_200cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
curve1 = interp_mean_d1_200cell' + interp_std_d1_200cell';
curve2 = interp_mean_d1_200cell' - interp_std_d1_200cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');
xlim([-1,5]);
yticks([0:0.2:0.8]);
ylim([0,0.8]);
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
legend([r3,b3],{'Non-persistent','Persistent'});
legend('boxoff');
title('Decoding with 200 Cells');

figure;
plot([50 100 200],first_point_delay(1,:), 'r-o','MarkerFaceColor','r','LineWidth',3)
hold on 
plot([50 100 200],first_point_delay(2,:), 'b-o','MarkerFaceColor','b','LineWidth',3)
xlabel('# Neurons')
ylabel('Accuracy')
box off
set(gca,'linewidth',1.5)
set(gca,'FontSize',16,'Fontweight','bold')
ylim([0.1,0.6]);

%%%%%%%%%%%%%%%%%%%%%%%%%overlay all%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
r1=plot(xq,interp_mean_d0_50cell,':r','LineWidth',3);
b1=plot(xq,interp_mean_d1_50cell,':b','LineWidth',3);
curve1 = interp_mean_d0_50cell' + interp_std_d0_50cell';
curve2 = interp_mean_d0_50cell' - interp_std_d0_50cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
curve1 = interp_mean_d1_50cell' + interp_std_d1_50cell';
curve2 = interp_mean_d1_50cell' - interp_std_d1_50cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');

r2=plot(xq,interp_mean_d0_100cell,'--r','LineWidth',3);
b2=plot(xq,interp_mean_d1_100cell,'--b','LineWidth',3);
curve1 = interp_mean_d0_100cell' + interp_std_d0_100cell';
curve2 = interp_mean_d0_100cell' - interp_std_d0_100cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
curve1 = interp_mean_d1_100cell' + interp_std_d1_100cell';
curve2 = interp_mean_d1_100cell' - interp_std_d1_100cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');

r3=plot(xq,interp_mean_d0_200cell,'-r','LineWidth',3);
b3=plot(xq,interp_mean_d1_200cell,'-b','LineWidth',3);
curve1 = interp_mean_d0_200cell' + interp_std_d0_200cell';
curve2 = interp_mean_d0_200cell' - interp_std_d0_200cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
curve1 = interp_mean_d1_200cell' + interp_std_d1_200cell';
curve2 = interp_mean_d1_200cell' - interp_std_d1_200cell';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');

xlim([-1,5]);
yticks([0:0.1:0.8]);
ylim([0,0.8]);
yline(0.125,'k','LineWidth',3);
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
%legend([r1,b1,r2,b2,r3,b3],{'50 cells Non-persistent','50 cells Persistent','100 cells Non-persistent','100 cells Persistent','200 cells Non-persistent','200 cells Persistent'});
%legend('boxoff');

figure;
hold on;
r1=plot(xq,interp_mean_d0_50cell_shuffle,':r','LineWidth',3);
b1=plot(xq,interp_mean_d1_50cell_shuffle,':b','LineWidth',3);
curve1 = interp_mean_d0_50cell_shuffle' + interp_std_d0_50cell_shuffle';
curve2 = interp_mean_d0_50cell_shuffle' - interp_std_d0_50cell_shuffle';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
curve1 = interp_mean_d1_50cell_shuffle' + interp_std_d1_50cell_shuffle';
curve2 = interp_mean_d1_50cell_shuffle' - interp_std_d1_50cell_shuffle';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');

r2=plot(xq,interp_mean_d0_100cell_shuffle,'--r','LineWidth',3);
b2=plot(xq,interp_mean_d1_100cell_shuffle,'--b','LineWidth',3);
curve1 = interp_mean_d0_100cell_shuffle' + interp_std_d0_100cell_shuffle';
curve2 = interp_mean_d0_100cell_shuffle' - interp_std_d0_100cell_shuffle';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
curve1 = interp_mean_d1_100cell_shuffle' + interp_std_d1_100cell_shuffle';
curve2 = interp_mean_d1_100cell_shuffle' - interp_std_d1_100cell_shuffle';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');

r3=plot(xq,interp_mean_d0_200cell_shuffle,'-r','LineWidth',3);
b3=plot(xq,interp_mean_d1_200cell_shuffle,'-b','LineWidth',3);
curve1 = interp_mean_d0_200cell_shuffle' + interp_std_d0_200cell_shuffle';
curve2 = interp_mean_d0_200cell_shuffle' - interp_std_d0_200cell_shuffle';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
curve1 = interp_mean_d1_200cell_shuffle' + interp_std_d1_200cell_shuffle';
curve2 = interp_mean_d1_200cell_shuffle' - interp_std_d1_200cell_shuffle';
x2 = [xq, fliplr(xq)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');

xlim([-1,5]);
yticks([0:0.1:0.8]);
ylim([0,0.8]);
yline(0.125,'k','LineWidth',3);
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
legend([r1,b1,r2,b2,r3,b3],{'50 cells Non-persistent','50 cells Persistent','100 cells Non-persistent','100 cells Persistent','200 cells Non-persistent','200 cells Persistent'});
legend('boxoff');
