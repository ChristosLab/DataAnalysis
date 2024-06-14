load('allcell_post.mat');  %analyse pre or post training
load('all_spatial_data.mat');
find_index=allcell_info{10};

temp_index=find_index(unique([allcell_info{3};allcell_info{4};allcell_info{5};allcell_info{6};allcell_info{7};allcell_info{8};allcell_info{9}]));

for r=1:10
    disp(r);
%if cell_group==1
%temp_vec=randperm(length(LS_index));
%temp_index=LS_index(temp_vec(1:num_select));
%else
%temp_vec=randperm(length(NMS_index));
%temp_index=NMS_index(temp_vec(1:num_select));
%end
remove_index=find(temp_index==1523);
temp_index(remove_index)=[];
remove_index=find(temp_index==1524);
temp_index(remove_index)=[];
remove_index=find(temp_index==255);
temp_index(remove_index)=[];
rng shuffle;
select_index=randperm(length(temp_index));
select_index=select_index(1:200);
[y,population_response]=numbi_normalize_svm_decoder(temp_index(select_index),all_spatial_data);
all_results(r,:)=y;
end

save postspatial_nonorm_numbi.mat all_results
