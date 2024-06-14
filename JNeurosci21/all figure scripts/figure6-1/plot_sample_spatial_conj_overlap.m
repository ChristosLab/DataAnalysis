%map cell category from passive to active task
load spatial_sample_f.mat
spatial_sample_f=f_sample;
load conj_sample_f.mat
conj_sample_f=f_sample_conj;
load find_index.mat
load spatial_inclusion_index.mat
num_cells=length(data_inclusion_index);
spatial_list=zeros(1,num_cells);
conj_list=zeros(1,num_cells);
load spatial_sample_info.mat;
spatial_info=category_sample_info;
spatial_CS=[spatial_info{5};spatial_info{6}];
[c,ia,ib]=intersect(data_inclusion_index,spatial_CS);
spatial_list(ia)=1;
spatial_NMS=[spatial_info{1};spatial_info{2};spatial_info{3};spatial_info{4}];
[c,ia,ib]=intersect(data_inclusion_index,spatial_NMS);
spatial_list(ia)=3;
spatial_LMS=[spatial_info{7}];
[c,ia,ib]=intersect(data_inclusion_index,spatial_LMS);
spatial_list(ia)=2;
all_spatial=[spatial_CS;spatial_NMS;spatial_LMS];
load conj_sample_info.mat;
conj_info=category_sample_info;
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

spatial_NMS_overlap=intersect(spatial_NMS,conj_NMS);
spatial_NMS_nonoverlap=setdiff(spatial_NMS,conj_NMS);

conj_NMS_overlap=intersect(conj_NMS,spatial_NMS);
conj_NMS_nonoverlap=setdiff(conj_NMS,spatial_NMS);

figure;
subplot(1,2,1);
hold on;
title('spatial cells');
scatter(ones(1,length(spatial_NMS_overlap)),spatial_sample_f(find_index(spatial_NMS_overlap),3));
scatter(2*ones(1,length(spatial_NMS_nonoverlap)),spatial_sample_f(find_index(spatial_NMS_nonoverlap),3));
subplot(1,2,2);
hold on;
title('conj cells');
scatter(ones(1,length(conj_NMS_overlap)),conj_sample_f(find_index(conj_NMS_overlap),3));
scatter(2*ones(1,length(conj_NMS_nonoverlap)),conj_sample_f(find_index(conj_NMS_nonoverlap),3));
%{
figure;
subplot(1,2,1)
bar([length(spatial_CS_overlap),length(spatial_CS)-length(spatial_CS_overlap);length(spatial_LMS_overlap),length(spatial_LMS)-length(spatial_LMS_overlap);length(spatial_NMS_overlap),length(spatial_NMS)-length(spatial_NMS_overlap)],'Stacked');
xticklabels({'CS','LMS','NMS'});
title('spatial');
ylim([0,120]);
subplot(1,2,2)
bar([length(conj_CS_overlap),length(conj_CS)-length(conj_CS_overlap);length(conj_LMS_overlap),length(conj_LMS)-length(conj_LMS_overlap);length(conj_NMS_overlap),length(conj_NMS)-length(conj_NMS_overlap)],'Stacked');
title('Conj');
ylim([0,120]);
xticklabels({'CS','LMS','NMS'});
all_overlap=intersect(all_conj,all_spatial);
%}
count1=length(intersect(find(spatial_list==3),find(conj_list==3)));
count2=length(intersect(find(spatial_list==3),find(conj_list~=3)));
count3=length(intersect(find(spatial_list~=3),find(conj_list==3)));
count4=length(intersect(find(spatial_list~=3),find(conj_list~=3)));
save all_overlap.mat all_overlap
save conj_list.mat conj_list
save spatial_list.mat spatial_list