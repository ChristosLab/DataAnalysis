load('allcell_pre.mat');  %analyse pre or post training
load('all_spatial_data.mat');
%cell_from_period=1; %1 from sample period, 2 from sampledelay period, 3 both period
find_index=allcell_info{10};
%if cell_from_period==1
%NMS_index=find_index(unique([allcell_info{3};allcell_info{4};allcell_info{5};allcell_info{6}]));
%LS_index=find_index(unique([allcell_info{7};allcell_info{8};allcell_info{9}]));
%num_select=min([length(NMS_index),length(LS_index)]);
%else
%NMS_index=find_index(unique([allcell_info{11};allcell_info{12};allcell_info{13};allcell_info{14}]));
%LS_index=find_index(unique([allcell_info{15};allcell_info{16};allcell_info{17}]));
%num_select=min([length(NMS_index),length(LS_index)]);
%end
%cell_group=2;  %1 represent LS, 2 represent NMS, 3 both group
temp_index=find_index(unique([allcell_info{3};allcell_info{4};allcell_info{5};allcell_info{6};allcell_info{7};allcell_info{8};allcell_info{9}]));
mode=3; %mode1 decode cue, mode2 decode sample, mode3 decode match/nonmatch 
remove_mode=1; %1 remove match pure, 2 remove location pure
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
%[x,y,population_response]=removepure_normalize_svm_decoder(0.4,0.1,temp_index(select_index),mode,remove_mode,all_spatial_data);
[x,y,population_response]=removepure_normalize_svm_decoder2(0.4,0.1,temp_index(select_index),mode,remove_mode,all_spatial_data);
all_results(r,:)=y;
end
saccuracy_mat=mean(all_results);
%save post_223_all.mat all_results
%save post_223_mean.mat accuracy_mat;
decode_std=std(all_results,0,1);
confidence_95=1.96*decode_std/sqrt(10);
x=linspace(-1,5,57);
xq=linspace(-1,5,400);
interp_mean=smooth(interp1(x,accuracy_mat,xq),10);
interp_95=smooth(interp1(x,confidence_95,xq),10);
figure;
plot(xq,interp_mean);
save post_all_match_orimean.mat interp_mean
save post_all_match_ori95.mat confidence_95
%save post_223_orimean.mat interp_mean
%save post_223_ori95.mat confidence_95