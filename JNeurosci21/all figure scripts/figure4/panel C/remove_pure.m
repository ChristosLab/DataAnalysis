function [new_mat]=remove_pure(old_mat, group_vec)
%old_mat is  num_condition cell, each cell contains spike count for trials within conditions, group_vet is condition group x column
%number of contions in each group in olad_mat matrix
%this version process data with 1 time point
%old_mat={[1.1,1.2,1.3,1.4],[2.2,2.3,2.4,2.5,2.6],[3.1,3.2,3.3],[4.1,4.2,4.3,4.4],[2.2,2.3,2.4,2.5,2.6],[3.1,3.2,3.3]};
%group_vec=[1,2;3,4;5,6];
num_groups=size(group_vec,1);
group_ids=[1:num_groups];
%condition_count=0;
for i=1:num_groups
    other_groups=setdiff(group_ids,i);%track current group under processing
    for j=1:size(group_vec,2) %loop through conditions in each group
 %       condition_count=condition_count+1;
        temp_targtrials=old_mat{group_vec(i,j)};%get spike time in selected condition
        targ_numtrial=length(temp_targtrials); %trial number in current condition
        group_select=[];
        for p=1:length(other_groups) %loop through other groups
            temp_conditions=group_vec(other_groups(p),:);
            source_pool=[];
            for q=1:length(temp_conditions)
                source_pool=[source_pool,old_mat{temp_conditions(q)}]; %pool all trials in other conditions in a big pool
            end
            group_select(p,:)=source_pool(randperm(length(source_pool),targ_numtrial));
        end
        new_mat{group_vec(i,j)}=mean([temp_targtrials;group_select]);
    end
end
end