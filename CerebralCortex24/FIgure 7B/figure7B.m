%delay interval for svm 200 neurons, 10 k fold, edges need to be 0.75 and
%1.75 for svm decoder
clear all 
close all

load('test_data.mat');  %analyse pre or post training
load('test_info.mat'); 
mode=1; %mode1 decode cue, mode2 decode sample, mode3 decode match/nonmatch 

%   standard SVM, n = 200, k = 10

for n = 1:100 %# of values we want to generate 
        disp(n);
    for r=1:3

          
        %D1=0
        temp_index = find(cell2mat(cellfun(@isnan,all_feature_info(:,2), 'UniformOutput', false))); %D1=NaN
        temp_index=temp_index(randperm(length(temp_index),180)); 
        [x,y,crv_cell,population_response]=normalize_SVM_decoder_delay(0.4,0.1,temp_index,mode);
        all_results0(r)=y;

        %d1=1
        temp_index = find(cell2mat(all_feature_info(:,2))==1);%D1=1
        temp_index=temp_index(randperm(length(temp_index),180)); 
        [x,y,crv_cell,population_response]=normalize_SVM_decoder_delay(0.4,0.1,temp_index,mode);
        all_results1(r)=y;

    end
    %d1=0
    accuracy_mat=mean(all_results0);
    delay_0(n) = accuracy_mat;
    

    %d1=1
    accuracy_mat=mean(all_results1);
    delay_1(n) = accuracy_mat;
   
end
 
save delay_d1=0_figure7 delay_0
save delay_d1=1_figure7 delay_1

%load('delay_d1=0_figure7.mat');
%load('delay_d1=1_figure7.mat');
figure;
hold on
b1=bar(categorical({'d1=0'}),mean(delay_0),0.2, 'red','EdgeColor','none');
b1.FaceAlpha=0.3;
b2=bar(categorical({'d1=1'}), mean(delay_1),0.2, 'blue','EdgeColor','none');
b2.FaceAlpha=0.3;
scatter(categorical({'d1=0'}), delay_0,50,"red")
scatter(categorical({'d1=1'}), delay_1,50, "blue")
ylabel('accuracy')
box off
set(gca,'linewidth',1)
set(gca,'FontSize',16,'Fontweight','bold')
hold off




