%spatial dimensionality err vs correct, location & match status trials
%peudopopulation was constructed once
load('all_spatial_info.mat');
load('all_spatial_data.mat');
load('all_spatialerr_info.mat');
load('all_spatialerr_data.mat');
group=all_spatial_info(:,4);  
IndexC = strfind(group,'POST');
find_index = find(not(cellfun('isempty',IndexC)));
for i=1:length(find_index)
    cell_correct={all_spatial_data{find_index(i),:}};
    cell_err={all_spatialerr_data{i,:}};
    for j=1:8
        temp_data=cell_correct{j};
        correct_count(i,(j-1)*2+1)=length(find([temp_data.IsMatch]==1));
        correct_count(i,(j-1)*2+2)=length(find([temp_data.IsMatch]==0));
        temp_data=cell_err{j};
        if length(temp_data)<1
           err_count(i,(j-1)*2+1)=0;
           err_count(i,(j-1)*2+2)=0;
        else
           if ~isfield(temp_data,'IsMatch')
              err_count(i,(j-1)*2+1)=0;
              err_count(i,(j-1)*2+2)=0;
           else
              err_count(i,(j-1)*2+1)=length(find([temp_data.IsMatch]==1));
              err_count(i,(j-1)*2+2)=length(find([temp_data.IsMatch]==0));
           end
        end
    end
end

match_select_vec=[1,3,5,7,9,11,13,15];%only matching trials analysed
num_condition_thred=[2:8];
num_trials_thred=[4:6];
for i=1:length(num_condition_thred)
    for j=1:length(num_trials_thred)
        pass_trial_requir=zeros(size(err_count,1),8);
        temp_condition_thred=num_condition_thred(i);
        temp_trials_thred=num_trials_thred(j);
        [temp_row,temp_col]=find(err_count(:,match_select_vec)>=temp_trials_thred & correct_count(:,match_select_vec)>=temp_trials_thred);
        for l=1:length(temp_row)
        pass_trial_requir(temp_row(l),temp_col(l))=1;
        end
        v = 1:1:8;
        C = nchoosek(v,temp_condition_thred);
        num_combination=size(C,1);
        cell_pass_count=[];
        for p=1:num_combination
           cell_pass_count(p)=sum(prod(pass_trial_requir(:,C(p,:)),2));
           cell_pass_index(p,:)=prod(pass_trial_requir(:,C(p,:)),2);
        end
        [pass_all_count(i,j),maximum_index]=max(cell_pass_count);
        select_condition{i,j}=C(maximum_index,:);
        select_cell{i,j}=find(cell_pass_index(maximum_index,:)==1);
    end
end

feed_matrix_correct=[];
feed_matrix_err=[];
select_thred=[2,1];
final_cells=select_cell{select_thred};
final_conditions=select_condition{select_thred};
mode=2; %1 for sample period 2 for sampledelay epriod
for i=1:length(final_cells)
    cell_correct={all_spatial_data{find_index(final_cells(i)),:}};
    cell_err={all_spatialerr_data{final_cells(i),:}};
    for j=1:length(final_conditions)
        class_correct=cell_correct{final_conditions(j)};
        %find match trials
        temp_index=find([class_correct.IsMatch]==1);

        if mode==1 %extract sample or sampledelay period
           temp_count=[class_correct.samplerate];
        else
           temp_count=[class_correct.sampledelay];
        end
        temp_count=temp_count(temp_index);
        randselect=randperm(length(temp_index));
        temp_count=temp_count(randselect(1:num_trials_thred(select_thred(2))));
        feed_matrix_correct(i,j,:)=temp_count;
        
        class_err=cell_err{final_conditions(j)};
        temp_index=find([class_err.IsMatch]==1);

        if mode==1 %extract sample or sampledelay period
           temp_count=[class_err.samplerate];
        else
           temp_count=[class_err.sampledelay];
        end
        temp_count=temp_count(temp_index);
        randselect=randperm(length(temp_index));
        temp_count=temp_count(randselect(1:num_trials_thred(select_thred(2))));
        feed_matrix_err(i,j,:)=temp_count;
        
    end
end

correct_winning_models=[];

for i=1:50
temp_index=randperm(size(feed_matrix_correct,1));
test_response=feed_matrix_correct(temp_index(1:50),:,:); %random select 100 cells
wholebrain_all{1}=test_response;
[subject_IDs, test_runs, winning_models, test_correlations]=functional_dimensionality_spike(wholebrain_all, ...
    ones(size(test_response,1),1),1,1);
correct_winning_models=[correct_winning_models,mean(winning_models)];
end
mean_correct_model = mean(correct_winning_models);
std_correct_model=1.96*std(correct_winning_models)/sqrt(50);

err_winning_models=[];

for i=1:50
temp_index=randperm(size(feed_matrix_err,1));
test_response=feed_matrix_err(temp_index(1:50),:,:); %random select 100 cells
wholebrain_all{1}=test_response;
[subject_IDs, test_runs, winning_models, test_correlations]=functional_dimensionality_spike(wholebrain_all, ...
    ones(size(test_response,1),1),1,1);
err_winning_models=[err_winning_models,mean(winning_models)];
end
mean_err_model = mean(err_winning_models);
std_err_model=1.96*std(err_winning_models)/sqrt(50);