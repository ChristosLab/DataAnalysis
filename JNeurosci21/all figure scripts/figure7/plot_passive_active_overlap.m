%map cell category from passive to active task
load passive_sample_f.mat
passive_f=f_sample;
load active_sample_f.mat
active_f=f_sample;
load active_trial96.mat
num_cells=length(index96);
passive_list=zeros(1,num_cells);
active_list=zeros(1,num_cells);
load passive_info.mat;
passive_info=category_info;
passive_CS=[passive_info{5};passive_info{6}];
passive_list(passive_CS)=1;
passive_NMS=[passive_info{1};passive_info{2};passive_info{3};passive_info{4}];
passive_list(passive_NMS)=3;
passive_LMS=[passive_info{7}];
passive_list(passive_LMS)=2;
all_passive=[passive_CS;passive_NMS;passive_LMS];
load active_info.mat;
active_info=category_info;
active_CS=[active_info{5};active_info{6}];
active_list(active_CS)=1;
active_NMS=[active_info{1};active_info{2};active_info{3};active_info{4}];
active_list(active_NMS)=3;
active_LMS=[active_info{7}];
active_list(active_LMS)=2;
all_active=[active_CS;active_NMS;active_LMS];

passive_NMS_overlap=intersect(passive_NMS,active_NMS);
passive_NMS_nonoverlap=setdiff(passive_NMS,active_NMS);

active_NMS_overlap=intersect(active_NMS,passive_NMS);
active_NMS_nonoverlap=setdiff(active_NMS,passive_NMS);

figure;
subplot(1,2,1);
hold on;
title('passive cells');
scatter(ones(1,length(passive_NMS_overlap)),passive_f(passive_NMS_overlap,3));
scatter(2*ones(1,length(passive_NMS_nonoverlap)),passive_f(passive_NMS_nonoverlap,3));
subplot(1,2,2);
hold on;
title('active cells');
scatter(ones(1,length(active_NMS_overlap)),active_f(active_NMS_overlap,3));
scatter(2*ones(1,length(active_NMS_nonoverlap)),active_f(active_NMS_nonoverlap,3));
%{
figure;
subplot(1,2,1)
bar([length(passive_CS_overlap),length(passive_CS)-length(passive_CS_overlap);length(passive_LMS_overlap),length(passive_LMS)-length(passive_LMS_overlap);length(passive_NMS_overlap),length(passive_NMS)-length(passive_NMS_overlap)],'Stacked');
xticklabels({'CS','LMS','NMS'});
title('Passive');
ylim([0,40]);
subplot(1,2,2)
bar([length(active_CS_overlap),length(active_CS)-length(active_CS_overlap);length(active_LMS_overlap),length(active_LMS)-length(active_LMS_overlap);length(active_NMS_overlap),length(active_NMS)-length(active_NMS_overlap)],'Stacked');
title('Active');
ylim([0,40]);
xticklabels({'CS','LMS','NMS'});
all_overlap=intersect(all_active,all_passive);
%}
count1=length(intersect(find(passive_list==3),find(active_list==3)));
count2=length(intersect(find(passive_list==3),find(active_list~=3)));
count3=length(intersect(find(passive_list~=3),find(active_list==3)));
count4=length(intersect(find(passive_list~=3),find(active_list~=3)));

save all_overlap.mat all_overlap
save active_list.mat active_list
save passive_list.mat passive_list