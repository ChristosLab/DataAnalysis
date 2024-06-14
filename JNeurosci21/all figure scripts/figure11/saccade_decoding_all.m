%decode saccade correct trial only, location & match status trials
%peudopopulation was constructed multiple times
%{
load('all_spatial_info.mat');
load('all_spatial_data.mat');
load('all_spatial_sacdirect.mat');
group=all_spatial_info(:,4);  
IndexC = strfind(group,'POST');
find_index = find(not(cellfun('isempty',IndexC)));
normalize=1;
for i=1:length(find_index)
    cell_correct={all_spatial_data{find_index(i),:}};
    cell_saccade=[all_trialdirect{find_index(i)}];
    if ~ismember(i,[1292,1293])
    if length(cell_saccade)>0
    for j=1:8  %loop through cell and condition, count trial numbers
        temp_data=cell_correct{j};
        correct_trialnum=[temp_data.trialnum];
        correct_class_sacdirect=cell_saccade(correct_trialnum);
        correct_count_sac1(i,(j-1)*2+1)=length(find([temp_data.IsMatch]==1 & correct_class_sacdirect==1));
        correct_count_sac2(i,(j-1)*2+1)=length(find([temp_data.IsMatch]==1 & correct_class_sacdirect==2));
        correct_count_sac1(i,(j-1)*2+2)=length(find([temp_data.IsMatch]==0 & correct_class_sacdirect==1));
        correct_count_sac2(i,(j-1)*2+2)=length(find([temp_data.IsMatch]==0 & correct_class_sacdirect==2));
    end
    end
   end
end


num_condition_thred=[16];
num_trials_thred=[3];
for i=1:length(num_condition_thred) %identify cells that quanlify
    for j=1:length(num_trials_thred)
        
        pass_trial_requir1=zeros(size(correct_count_sac1));
        pass_trial_requir2=zeros(size(correct_count_sac2));
        temp_condition_thred=num_condition_thred(i);
        temp_trials_thred=num_trials_thred(j);
        [temp_row1,temp_col1]=find(correct_count_sac1>=temp_trials_thred);
        [temp_row2,temp_col2]=find(correct_count_sac2>=temp_trials_thred);
        for L=1:length(temp_row1)
        pass_trial_requir1(temp_row1(L),temp_col1(L))=1;
        end
        for L=1:length(temp_row2)
        pass_trial_requir2(temp_row2(L),temp_col2(L))=1;
        end
        pass_trial_require=pass_trial_requir1.*pass_trial_requir2;
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
%}
load('post_spatial_inclusion.mat');
correctcell_index=post_inclusion;
%use only certain category
load('choice_CS.mat');
temp1=temp;
load('choice_LMS.mat');
temp2=temp;
temp=unique([temp1;temp2]);
%load('choice_NMS.mat');
cell_pool=intersect(correctcell_index,temp);

for p=1:10
    disp(p);
    randsample=cell_pool(randperm(length(cell_pool),50));% randomly choose 50
   %[x,y]=spatial_saccade_decoder_all2(0.4,0.1,correctcell_index); %both versions of spatial_saccade_decoder_all work
   [x,y]=spatial_saccade_decoder_all2(0.4,0.1,randsample); %both versions of spatial_saccade_decoder_all work
   all_results(p,:)=y;
end
accuracy_mat=mean(all_results);
%save decode_sac_all.mat all_results
%save decode_sac_mean.mat accuracy_mat;
save decode_sac_linear50.mat all_results
save decode_sac_linearmean.mat accuracy_mat;
decode_std=std(all_results,0,1);
confidence_95=1.96*decode_std/sqrt(10);
x=linspace(-1,6,67);
xq=linspace(-1,6,400);
interp_mean=smooth(interp1(x,accuracy_mat,xq),10);
interp_95=smooth(interp1(x,confidence_95,xq),10);
figure;
plot(xq,interp_mean);
