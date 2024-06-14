load('all_feature_data.mat');
load('all_feature_info.mat');
load('passive_informative.mat');
passive_informative=find_informative;
load('active_informative.mat');
active_informative=find_informative;
load passive_info.mat; %dataset to use
passive_info=category_info;
load active_info.mat;
active_info=category_info;

stim_set=2;  %which stimuli set to use 

unique_sets=unique(all_feature_info(:,3));
feature10_group=[2,3,4,5,6,7];
featurex_group=[1,9,11,14,19,22,24];
featurexa_group=[8,10,12,15,16,18,20,23,25,26];
featurexb_group=[13,17,21,27,28];

grouplookup{1}=[2,3,4,5,6,7,8,1];
grouplookup{2}=[2,3,8,5,6,7,4,1];
grouplookup{3}=[5,3,8,2,6,7,4,1];
grouplookup{4}=[2,3,4,5,6,7,8,1,10,9];
lookup_table=grouplookup{stim_set};
for i=1:size(all_feature_info,1)
    cell_set=all_feature_info{i,3};
    cell_setnum=find(strcmp(unique_sets,cell_set));
    if ismember(cell_setnum,feature10_group)
        cell_group_num(i)=4;
    elseif ismember(cell_setnum,featurex_group)
        cell_group_num(i)=1;
    elseif ismember(cell_setnum,featurexa_group)
        cell_group_num(i)=2;
    else
        cell_group_num(i)=3;
    end
end
cell_from_period=2; %1 from sample period, 2 from sampledelay period

%cells_to_useall=passive_informative;
%cells_to_useall=active_informative;
if cell_from_period==1   %decode sensory and sensory use correponding cells
pre_LS=unique([passive_info{5};passive_info{6};passive_info{7}]);
post_LS=unique([active_info{5};active_info{6};active_info{7}]);
pre_NMS=unique([passive_info{1};passive_info{2};passive_info{3};passive_info{4}]);
post_NMS=unique([active_info{1};active_info{2};active_info{3};active_info{4}]);
%pre_LS_stimset=intersect(pre_LS,find(cell_group_num==stim_set));
%post_LS_stimset=intersect(post_LS,find(cell_group_num==stim_set));
%pre_NMS_stimset=intersect(pre_NMS,find(cell_group_num==stim_set));
%post_NMS_stimset=intersect(post_NMS,find(cell_group_num==stim_set));
pre_LS_stimset=pre_LS;
post_LS_stimset=post_LS;
pre_NMS_stimset=pre_NMS;
post_NMS_stimset=post_NMS;
pre_all=unique([pre_LS_stimset;pre_NMS_stimset]);
post_all=unique([post_LS_stimset;post_NMS_stimset]);
numpick=min([length(pre_all),length(post_all)]);
mode=2;
else
pre_LS=unique([passive_info{12};passive_info{13};passive_info{14}]);
post_LS=unique([active_info{12};active_info{13};active_info{14}]);
pre_NMS=unique([passive_info{8};passive_info{9};passive_info{10};passive_info{11}]);
post_NMS=unique([active_info{8};active_info{9};active_info{10};active_info{11}]);
%pre_LS_stimset=intersect(pre_LS,find(cell_group_num==stim_set));
%post_LS_stimset=intersect(post_LS,find(cell_group_num==stim_set));
%pre_NMS_stimset=intersect(pre_NMS,find(cell_group_num==stim_set));
%post_NMS_stimset=intersect(post_NMS,find(cell_group_num==stim_set));
pre_LS_stimset=pre_LS;
post_LS_stimset=post_LS;
pre_NMS_stimset=pre_NMS;
post_NMS_stimset=post_NMS;
pre_all=unique([pre_LS_stimset;pre_NMS_stimset]);
post_all=unique([post_LS_stimset;post_NMS_stimset]);
numpick=min([length(pre_all),length(post_all)]);
mode=3;
end

for r=1:10  %repeat time
    disp(r);
%numpick=150;
%temp_vec=randperm(length(pre_all));
%cells_to_use_stimset=pre_all(temp_vec(1:numpick)); 
%remove_index=find(cells_to_use_stimset==196);
%cells_to_use_stimset(remove_index)=[];
%remove_index=find(cells_to_use_stimset==167);
%cells_to_use_stimset(remove_index)=[];
%remove_index=find(cells_to_use_stimset==194);
%cells_to_use_stimset(remove_index)=[];
%remove_index=find(cells_to_use_stimset==195);
%cells_to_use_stimset(remove_index)=[];

numpick=150;
temp_vec=randperm(length(post_all));
cells_to_use_stimset=post_all(temp_vec(1:numpick)); 

%[x,y,crv_cell,population_response]=normalize_SVM_decoder3(0.4,0.1,cells_to_use_stimset,mode,stim_set); 
[x,y,crv_cell,population_response]=normalize_SVM_decoder3(0.4,0.1,cells_to_use_stimset+1553,mode,stim_set); 
store_mat(r,:)=y;
end
disp(size(cells_to_use_stimset));
plot_mat=mean(store_mat,1);
save pre_decodematch_all.mat store_mat
save pre_decodematch_mean.mat plot_mat;
decode_std=std(store_mat,0,1);
confidence_95=1.96*decode_std/sqrt(10);
x=linspace(-1,5,57);
xq=linspace(-1,5,400);
interp_mean=smooth(interp1(x,plot_mat,xq),10);
plot(xq,interp_mean);
save pre_decodematch_orimean.mat interp_mean
save pre_decodematch_ori95.mat confidence_95