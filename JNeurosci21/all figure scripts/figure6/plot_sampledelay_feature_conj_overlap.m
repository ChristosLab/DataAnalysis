%map cell category from passive to active task
load feature_sampledelay_f.mat
feature_sampledelay_f=f_sampledelay;
load conj_sampledelay_f.mat
conj_sampledelay_f=f_sampledelay_conj;
load find_index.mat
load feature_inclusion_index.mat
num_cells=length(data_inclusion_index);
feature_list=zeros(1,num_cells);
conj_list=zeros(1,num_cells);
load feature_sampledelay_info.mat;
feature_info=category_sampledelay_info;
feature_CS=[feature_info{5};feature_info{6}];
[c,ia,ib]=intersect(data_inclusion_index,feature_CS);
feature_list(ia)=1;
feature_NMS=[feature_info{1};feature_info{2};feature_info{3};feature_info{4}];
[c,ia,ib]=intersect(data_inclusion_index,feature_NMS);
feature_list(ia)=3;
feature_LMS=[feature_info{7}];
[c,ia,ib]=intersect(data_inclusion_index,feature_LMS);
feature_list(ia)=2;
all_feature=[feature_CS;feature_NMS;feature_LMS];
load conj_sampledelay_info.mat;
conj_info=category_sampledelay_info;
conj_CS=[conj_info{5};conj_info{6}];
[c,ia,ib]=intersect(data_inclusion_index,conj_CS);
conj_list(ia)=1;
conj_NMS=[conj_info{1};conj_info{2};conj_info{3};conj_info{4}];
[c,ia,ib]=intersect(data_inclusion_index,conj_NMS);
conj_list(ia)=3;
conj_LMS=[conj_info{7}];
[c,ia,ib]=intersect(data_inclusion_index,conj_LMS);
conj_list(ia)=2;
all_conj=[conj_CS;conj_NMS;conj_LMS];


feature_NMS_overlap=intersect(feature_NMS,conj_NMS);
feature_NMS_nonoverlap=setdiff(feature_NMS,conj_NMS);


conj_NMS_overlap=intersect(conj_NMS,feature_NMS);
conj_NMS_nonoverlap=setdiff(conj_NMS,feature_NMS);


figure;
subplot(1,2,1);
hold on;
title('feature cells');
scatter(ones(1,length(feature_NMS_overlap)),feature_sampledelay_f(find_index(feature_NMS_overlap),3));
scatter(2*ones(1,length(feature_NMS_nonoverlap)),feature_sampledelay_f(find_index(feature_NMS_nonoverlap),3));
subplot(1,2,2);
hold on;
title('conj cells');
scatter(ones(1,length(conj_NMS_overlap)),conj_sampledelay_f(find_index(conj_NMS_overlap),3));
scatter(2*ones(1,length(conj_NMS_nonoverlap)),conj_sampledelay_f(find_index(conj_NMS_nonoverlap),3));
%{
figure;
subplot(1,2,1)
bar([length(feature_CS_overlap),length(feature_CS)-length(feature_CS_overlap);length(feature_LMS_overlap),length(feature_LMS)-length(feature_LMS_overlap);length(feature_NMS_overlap),length(feature_NMS)-length(feature_NMS_overlap)],'Stacked');
xticklabels({'CS','LMS','NMS'});
title('feature');
ylim([0,120]);
subplot(1,2,2)
bar([length(conj_CS_overlap),length(conj_CS)-length(conj_CS_overlap);length(conj_LMS_overlap),length(conj_LMS)-length(conj_LMS_overlap);length(conj_NMS_overlap),length(conj_NMS)-length(conj_NMS_overlap)],'Stacked');
title('Conj');
ylim([0,120]);
xticklabels({'CS','LMS','NMS'});
all_overlap=intersect(all_conj,all_feature);
%}
count1=length(intersect(find(feature_list==3),find(conj_list==3)));
count2=length(intersect(find(feature_list==3),find(conj_list~=3)));
count3=length(intersect(find(feature_list~=3),find(conj_list==3)));
count4=length(intersect(find(feature_list~=3),find(conj_list~=3)));
%save all_overlap.mat all_overlap
%save conj_list.mat conj_list
%save feature_list.mat feature_list