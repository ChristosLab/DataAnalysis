clear all 
close all

load('feature_sig_data.mat');  %analyse pre or post training
load('feature_sig_info.mat'); 
load('feature_nonsig_data.mat'); 
mode=1; %mode1 decode cue, mode2 decode sample, mode3 decode match/nonmatch 

%   standard SVM, n = 200, k = 10
interp_mean0 = []; %sig-nonpersistent
interp_mean1 = []; %persistent
interp_mean2 = []; %all_nonpersistent
for n = 1:100

   % for r=1:3
        disp(n);
          
        %D1=0
        temp_index = find(cell2mat(cellfun(@isnan,all_feature_info(:,2), 'UniformOutput', false))); %D1=NaN
        temp_index=temp_index(randperm(length(temp_index),180)); 
        [x,y,crv_cell,population_response]=normalize_SVM_decoder_providedata(0.4,0.1,temp_index,all_feature_data,mode);
        all_results0=y;

        temp_index = find(cell2mat(all_feature_info(:,2))==1);%D1=1
        temp_index=temp_index(randperm(length(temp_index),180)); 
        [~,y,crv_cell,population_response]=normalize_SVM_decoder_providedata(0.4,0.1,temp_index,all_feature_data,mode);
        all_results1=y;
        
        %all non-persistent inclcuding D1==-1
        temp_index = find(cell2mat(all_feature_info(:,2))~=1);%D1=1
        temp_index=[temp_index;[size(all_feature_data,1)+1:size(all_feature_data,1)+size(nonsig_feature_data,1)]'];
        temp_index=temp_index(randperm(length(temp_index),180)); 
        [x,y,crv_cell,population_response]=normalize_SVM_decoder_providedata(0.4,0.1,temp_index,[all_feature_data;nonsig_feature_data],mode);
        all_results2=y;


   % end

    %d1=0
    accuracy_mat=all_results0;
    x=linspace(-1.1,5,length(y));%57 
    xq=linspace(-1.1,5,400);
    interp_mean0(n, :)=smooth(interp1(x,accuracy_mat,xq),10);

    %d1=1
    accuracy_mat=all_results1;
    x=linspace(-1.1,5,length(y));%57 
    xq=linspace(-1.1,5,400);
    interp_mean1(n, :)=smooth(interp1(x,accuracy_mat,xq),10);
    
    accuracy_mat=all_results2;
    x=linspace(-1.1,5,length(y));%57 
    xq=linspace(-1.1,5,400);
    interp_mean2(n, :)=smooth(interp1(x,accuracy_mat,xq),10);

   
end

save feature_sig_nonpersistent_figure2 interp_mean0
save feature_persistent_figure2 interp_mean1
save feature_all_nonpersistent_figure2 interp_mean2

load('feature_sig_nonpersistent_figure2.mat');
load('feature_persistent_figure2.mat');
load('feature_all_nonpersistent_figure2.mat');


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
plot(x, y, 'r', 'LineWidth', 2);
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
xlim([-1,5]);

%d1=1
y = mean(interp_mean1); 
x = xq;
std_dev = std(interp_mean1,0,1);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b', 'FaceAlpha',0.3, 'EdgeColor','none');
hold on;
plot(x, y, 'b', 'LineWidth', 2);
xlabel('Time (s)')
ylabel('Accuracy')
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
xlim([-1,5]);
%{
y = mean(interp_mean2); 
x = xq;
std_dev = std(interp_mean2,0,1);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'g', 'FaceAlpha',0.3, 'EdgeColor','none');

plot(x, y, 'g', 'LineWidth', 2);
xlabel('Time (s)')
ylabel('Accuracy')
%}

hold off
