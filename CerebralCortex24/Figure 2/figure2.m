clear all 
close all

load('spatial_sig_data.mat');  %analyse pre or post training
load('spatial_sig_info.mat'); 
load('spatial_nonsig_data.mat'); 
mode=1; %mode1 decode cue, mode2 decode sample, mode3 decode match/nonmatch 

%   standard SVM, n = 200, k = 10
interp_mean0 = []; %sig-nonpersistent
interp_mean1 = []; %persistent
interp_mean2 = []; %all_nonpersistent
interp_mean3 = []; %all_inhib
interp_mean4 = []; %persistent 50 cells
for n = 1:100

   % for r=1:3
        disp(n);
        
        %D1=0
        temp_index = find(cell2mat(cellfun(@isnan,all_spatial_info(:,2), 'UniformOutput', false))); %D1=NaN
        temp_index=temp_index(randperm(length(temp_index),200)); 
        [x,y,crv_cell,population_response]=normalize_SVM_decoder_providedata(0.4,0.1,temp_index,all_spatial_data,mode);
        all_results0=y;

        temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
        temp_index=temp_index(randperm(length(temp_index),200)); 
        [~,y,crv_cell,population_response]=normalize_SVM_decoder_providedata(0.4,0.1,temp_index,all_spatial_data,mode);
        all_results1=y;
        
        %all non-persistent
        temp_index = find(cell2mat(all_spatial_info(:,2))~=1);%D1=-1
        temp_index=[temp_index;[size(all_spatial_data,1)+1:size(all_spatial_data,1)+size(nonsig_spatial_data,1)]'];
        temp_index=temp_index(randperm(length(temp_index),200)); 
        [x,y,crv_cell,population_response]=normalize_SVM_decoder_providedata(0.4,0.1,temp_index,[all_spatial_data;nonsig_spatial_data],mode);
        all_results2=y;

        %all inhib
        temp_index = find(cell2mat(all_spatial_info(:,2))==-1);%D1=-1
        temp_index=temp_index(randperm(length(temp_index),50)); 
        [~,y,crv_cell,population_response]=normalize_SVM_decoder_providedata(0.4,0.1,temp_index,all_spatial_data,mode);
        all_results3=y;
        
        temp_index = find(cell2mat(all_spatial_info(:,2))==1);%D1=1
        temp_index=temp_index(randperm(length(temp_index),50)); 
        [~,y,crv_cell,population_response]=normalize_SVM_decoder_providedata(0.4,0.1,temp_index,all_spatial_data,mode);
        all_results4=y;
   % end

    %d1=0
    accuracy_mat=all_results0;
    x=linspace(-1.1,5.4,length(x));%57 
    xq=linspace(-1.1,5.4,400);
    interp_mean0(n, :)=smooth(interp1(x,accuracy_mat,xq),10);

    %d1=1
    accuracy_mat=all_results1;
    x=linspace(-1.1,5.4,length(x));%57 
    xq=linspace(-1.1,5.4,400);
    interp_mean1(n, :)=smooth(interp1(x,accuracy_mat,xq),10);
    
    accuracy_mat=all_results2;
    x=linspace(-1.1,5.4,length(x));%57 
    xq=linspace(-1.1,5.4,400);
    interp_mean2(n, :)=smooth(interp1(x,accuracy_mat,xq),10);

    accuracy_mat=all_results3;
    x=linspace(-1.1,5.4,length(x));%57 
    xq=linspace(-1.1,5.4,400);
    interp_mean3(n, :)=smooth(interp1(x,accuracy_mat,xq),10);
    
    accuracy_mat=all_results4;
    x=linspace(-1.1,5.4,length(x));%57 
    xq=linspace(-1.1,5.4,400);
    interp_mean4(n, :)=smooth(interp1(x,accuracy_mat,xq),10);
end

save spatial_sig_nonpersistent_figure2 interp_mean0
save spatial_persistent_figure2 interp_mean1
save spatial_all_nonpersistent_figure2 interp_mean2
save spatial_all_inhib_figure2 interp_mean3
save spatial_persistent50_figure2 interp_mean4

load spatial_sig_nonpersistent_figure2 
load spatial_persistent_figure2 
load spatial_all_nonpersistent_figure2 
load spatial_all_inhib_figure2 
load spatial_persistent50_figure2 

%Figure 2b
figure;
xq=linspace(-1.1,5,400);
%d1=0
y = mean(interp_mean0); 
x = xq;
std_dev = std(interp_mean0,0,1);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
hold on;
l2=plot(x, y, 'r', 'LineWidth', 2);
%d1=1
y = mean(interp_mean1); 
x = xq;
std_dev = std(interp_mean1,0,1);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');
l1=plot(x, y, 'b', 'LineWidth', 2);
xlabel('Time (s)')
ylabel('Accuracy')
xlim([-1,5]);
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
legend( [l1;l2] , {'Persistent','Non-persistent'} );
legend('boxoff');
%Figure 2C
figure;
hold on;
y = mean(interp_mean1); 
x = xq;
std_dev = std(interp_mean1,0,1);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');
l1=plot(x, y, 'b', 'LineWidth', 2);

y = mean(interp_mean2); 
x = xq;
std_dev = std(interp_mean2,0,1);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
l2=plot(x, y, 'r', 'LineWidth', 2);
xlabel('Time (s)')
ylabel('Accuracy')
xlim([-1,5]);
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
legend( [l1;l2] , {'Persistent','Task responsive non-persistent'} );
legend('boxoff');
%Figure 2D
figure;
hold on;
y = mean(interp_mean4); 
x = xq;
std_dev = std(interp_mean1,0,1);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');
l1=plot(x, y, 'b', 'LineWidth', 2);

y = mean(interp_mean3); 
x = xq;
std_dev = std(interp_mean2,0,1);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.3, 'EdgeColor','none');
l2=plot(x, y, 'r', 'LineWidth', 2);
xlabel('Time (s)')
ylabel('Accuracy')
xlim([-1,5]);
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
legend( [l1;l2] , {'Persistent','Inhibited'} );
legend('boxoff');