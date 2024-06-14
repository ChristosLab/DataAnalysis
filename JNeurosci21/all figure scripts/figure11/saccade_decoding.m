%decode saccade err vs correct, location & match status trials
%peudopopulation was constructed multiple times
load('all_spatial_info.mat');
load('all_spatial_data.mat');
load('all_spatialerr_info.mat');
load('all_spatialerr_data.mat');
load('spatial_sacdirect.mat');
group=all_spatial_info(:,4);  
IndexC = strfind(group,'POST');
load('post_spatial_inclusion.mat');
find_index = post_inclusion;
normalize=1;
for i=1:length(find_index)
    cell_correct={all_spatial_data{find_index(i),:}};
    cell_err={all_spatialerr_data{find_index(i)-1735,:}};
    cell_saccade=[all_trialdirect{find_index(i)}];
    cell_correct_sacdirect=[];
    cell_err_sacdirect=[];
    cell_correct_ismatch=[];
    cell_err_ismatch=[];

  %  if ~ismember(i,[])
     for j=1:8  %loop cuecalss
        temp_correctdata=cell_correct{j};%get physioogy data
        correct_trialnum=[temp_correctdata.trialnum]; %correct trial index
        correct_class_sacdirect=cell_saccade(correct_trialnum); %correct trial sac direction
        cell_correct_ismatch=[cell_correct_ismatch,[temp_correctdata.IsMatch]];
        cell_correct_sacdirect=[cell_correct_sacdirect,correct_class_sacdirect];
        temp_errdata=cell_err{j};%get physioogy data
        if length(temp_errdata)>0
            if isfield(temp_errdata,'trialnum')
            err_trialnum=[temp_errdata.trialnum];
            else
            err_trialnum=[temp_errdata.Trial_Num];
            end
        err_class_sacdirect=cell_saccade(err_trialnum); %correct trial sac direction
        cell_err_ismatch=[cell_err_ismatch,[temp_errdata.IsMatch]];
        cell_err_sacdirect=[cell_err_sacdirect,err_class_sacdirect];
        else
        cell_err_ismatch=[];
        cell_err_sacdirect=[];
        end
     end
     for j=1:8
     correct_count_sac(i,j)=length(find(cell_correct_sacdirect==j & cell_correct_ismatch==1));
     correct_count_sac(i,j+8)=length(find(cell_correct_sacdirect==j & cell_correct_ismatch==0));
     err_count_sac(i,j)=length(find(cell_err_sacdirect==j & cell_err_ismatch==1));
     err_count_sac(i,j+8)=length(find(cell_err_sacdirect==j & cell_err_ismatch==0));
     end
 %end
end


num_condition_thred=[2:8];
num_trials_thred=[1:4];
for i=1:length(num_condition_thred) %identify cells that quanlify
    for j=1:length(num_trials_thred)       
        pass_trial_require=zeros(size(err_count_sac));
        temp_condition_thred=num_condition_thred(i);
        temp_trials_thred=num_trials_thred(j);
        [temp_row,temp_col]=find(err_count_sac>=temp_trials_thred & correct_count_sac>=temp_trials_thred);
        for L=1:length(temp_row)
        pass_trial_require(temp_row(L),temp_col(L))=1;
        end
        v = 1:1:16;
        C = nchoosek(v,temp_condition_thred);
        num_combination=size(C,1);
        cell_pass_count=[];
        for p=1:num_combination
           cell_pass_count(p)=sum(prod(pass_trial_require(:,C(p,:)),2));
           cell_pass_index(p,:)=prod(pass_trial_require(:,C(p,:)),2);
        end
        [pass_all_count(i,j),maximum_index]=max(cell_pass_count);
        select_condition{i,j}=C(maximum_index,:);
        select_cell{i,j}=find(cell_pass_index(maximum_index,:)==1);        
    end
end
cell_list=select_cell{1,2};
condition_list=select_condition{1,2};
%only use CS+LMS or NMS
%load('choice_CS.mat');
%temp1=temp;
%load('choice_LMS.mat');
%temp2=temp;
%temp=unique([temp1;temp2]);
load('choice_NMS.mat');
category_cell=intersect(find_index(cell_list),temp);

for p=1:10
    disp(p);
    category_cell=category_cell(randperm(length(category_cell),10));
   [x,y1,y2,yall]=spatial_saccade_decoder(0.4,0.1,category_cell,category_cell-1735,condition_list);
  % [x,y1,y2,yall]=spatial_saccade_decoder(0.4,0.1,find_index(cell_list),find_index(cell_list)-1735,condition_list);
   correct_results(p,:)=y1;
   err_results(p,:)=y2;
   allcell_results(p,:)=yall;
end

correct_accuracy_mat=mean(correct_results);
err_accuracy_mat=mean(err_results);
allcell_accuracy_mat=mean(allcell_results);
all_results=[correct_results;err_results;allcell_results];
accuracy_mat=[correct_accuracy_mat;err_accuracy_mat;allcell_accuracy_mat];
save sac_correctvserr_nonlinearall.mat all_results
save sac_correctvserr_nonlinearmean.mat accuracy_mat;
%decode_std=std(all_results,0,1);
%confidence_95=1.96*decode_std/sqrt(10);
x=linspace(-1,6,67);
xq=linspace(-1,6,400);
interp_mean_correct=smooth(interp1(x,correct_accuracy_mat,xq),10);
interp_mean_err=smooth(interp1(x,err_accuracy_mat,xq),10);

%interp_95=smooth(interp1(x,confidence_95,xq),10);
figure;
plot(xq,interp_mean_correct);
hold on;
plot(xq,interp_mean_err);
%{
%for plotting
results=all_results(11:20,:);
condition_mean=accuracy_mat(2,:);
decode_std=std(results,0,1);
confidence_95=1.96*decode_std/sqrt(10);
x=linspace(-1,6,67);
xq=linspace(-1,6,400);
interp_mean=smooth(interp1(x,condition_mean,xq),10);
interp_95=smooth(interp1(x,confidence_95,xq),10);

%plot correct vs err
load('sac_correctvserr_linearall.mat')
load('sac_correctvserr_linearmean.mat');
figure;
plot(accuracy_mat(1,:));
hold on;
plot(accuracy_mat(2,:));
%}