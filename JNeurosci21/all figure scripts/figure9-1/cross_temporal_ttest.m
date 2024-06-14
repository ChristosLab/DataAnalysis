load('post_delay_1_3.mat');
S1=all_results;
load('post_delay_2_3.mat');
S2=all_results;
for i=1:size(S1,2)
    for j=1:size(S2,3)
    [H,P]=ttest2(S1(:,i,j),S2(:,i,j));
    ori_matrix(i,j)=P;
    end
end
orithred_matrix=zeros(size(ori_matrix));
orithred_matrix(ori_matrix<0.05)=1;

rng('shuffle');
shuffle_pool=[S1;S2];
num_pertubation=10;
for p=1:num_pertubation
    temp_pertubation=randperm(20);
    shuffle_S1=shuffle_pool(temp_pertubation(1:10),:,:);
    shuffle_S2=shuffle_pool(temp_pertubation(11:20),:,:);
  for i=1:size(shuffle_S1,2)
    for j=1:size(shuffle_S2,3)
    [H,P]=ttest2(shuffle_S1(:,i,j),shuffle_S2(:,i,j));
    permute_matrix(i,j)=P;
    end
  end
permute_extreme=min(min(permute_matrix));  %get streme stats value for pixel based correction
pixel_extreme(p)=permute_extreme;

thred_matrix=zeros(size(permute_matrix));
thred_matrix(permute_matrix<0.05)=1;
[cluster_matrix,cluster_num]=bwlabeln(thred_matrix);  %label all clusters that is below p value of 0.05
biggest_cluster=0;
  for i=1:cluster_num                                   %get the size of biggest cluster in each interation
    cluster_size=length(find(cluster_matrix==i));
    if cluster_size>biggest_cluster
        biggest_cluster=cluster_size;
    end
  end
extreme_cluster_size(p)=biggest_cluster;
end

 [cluster_matrix,cluster_num]=bwlabeln(orithred_matrix);  %find original significant clusters 
 for i=1:cluster_num
     size_cluster=length(find(cluster_matrix==i)); %find size and location for cluster in original comparison
     location_cluster=find(cluster_matrix==i);
     if sum(extreme_cluster_size>size_cluster)/length(extreme_cluster_size)>0.05
         cluster_matrix(location_cluster)=0;
     end  %find the clusters in original comparison that is larger than chance
 end

 
accuracy_S1=squeeze(mean(S2));
[X,Y] = meshgrid(-0.8:0.1:4.8);
[Xq,Yq] = meshgrid(-0.8:0.02:4.8);
Vq = interp2(X,Y,accuracy_S1,Xq,Yq);
Vq=Vq-median(reshape(Vq,1,size(Vq,1)*size(Vq,2)));
Vq=Vq/max(max(Vq));
imagesc(flipud(Vq));
colormap parula;
caxis([0,1]);
hold on;
plot([40,40],[0,281],'k','LineWidth',2);
plot([65,65],[0,281],'k','LineWidth',2);
plot([140,140],[0,281],'k','LineWidth',2);
plot([165,165],[0,281],'k','LineWidth',2);
plot([0,281],[240,240],'k','LineWidth',2);
plot([0,281],[215,215],'k','LineWidth',2);
plot([0,281],[140,140],'k','LineWidth',2);
plot([0,281],[115,115],'k','LineWidth',2);
yticks([40 90 140 190 240]);
yticklabels({'4','3','2','1','0'});
a = get(gca,'YTickLabel'); 
set(gca,'YTickLabel',a,'fontsize',20);
xticks([40 90 140 190 240]);
xticklabels({'0','1','2','3','4'});
a = get(gca,'XTickLabel'); 
set(gca,'XTickLabel',a,'fontsize',20);

cluster_matrix(cluster_matrix>0)=1;
C=bwboundaries(cluster_matrix);
%{
 for k = 1:length(C)
     boundary = C{k};
     interp_boundx=5*(boundary(:,1)-1);
     interp_boundy=282-5*(boundary(:,2)-1);
     plot(interp_boundx, interp_boundy, 'g', 'LineWidth', 2);
 end
%}
binary_sig_stats=ori_matrix<min(pixel_extreme);  %because we permute for 20 times, 0.05 correspond to smaller than the smallest in the extreme value vector, if we permute 200 times 0.05 will correspond to smaller than the 10th smallest value in extreme stats vector
B = bwboundaries(binary_sig_stats);
for k = 1:length(B)
    boundary_size(k)=size(B{k},1);
end
[sort_size,I]=sort(boundary_size,'descend');
for k = 1:1
    boundary = B{I(k)};
    interp_boundx=5*(boundary(:,1)-1);
    interp_boundy=282-5*(boundary(:,2)-1);
    plot(interp_boundx, interp_boundy, 'r', 'LineWidth', 2)
end

figure;
accuracy_S1=squeeze(mean(S1));
[X,Y] = meshgrid(-0.8:0.1:4.8);
[Xq,Yq] = meshgrid(-0.8:0.02:4.8);
Vq = interp2(X,Y,accuracy_S1,Xq,Yq);
Vq=Vq-median(reshape(Vq,1,size(Vq,1)*size(Vq,2)));
Vq=Vq/max(max(Vq));
imagesc(flipud(Vq));
colormap parula;
caxis([0,1]);
hold on;
plot([40,40],[0,281],'k','LineWidth',2);
plot([65,65],[0,281],'k','LineWidth',2);
plot([140,140],[0,281],'k','LineWidth',2);
plot([165,165],[0,281],'k','LineWidth',2);
plot([0,281],[240,240],'k','LineWidth',2);
plot([0,281],[215,215],'k','LineWidth',2);
plot([0,281],[140,140],'k','LineWidth',2);
plot([0,281],[115,115],'k','LineWidth',2);
yticks([40 90 140 190 240]);
yticklabels({'4','3','2','1','0'});
a = get(gca,'YTickLabel'); 
set(gca,'YTickLabel',a,'fontsize',20);
xticks([40 90 140 190 240]);
xticklabels({'0','1','2','3','4'});
a = get(gca,'XTickLabel'); 
set(gca,'XTickLabel',a,'fontsize',20);

cluster_matrix(cluster_matrix>0)=1;
C=bwboundaries(cluster_matrix);
%{
 for k = 1:length(C)
     boundary = C{k};
     interp_boundx=5*(boundary(:,1)-1);
     interp_boundy=282-5*(boundary(:,2)-1);
     plot(interp_boundx, interp_boundy, 'g', 'LineWidth', 2);
 end
%}
binary_sig_stats=ori_matrix<min(pixel_extreme);  %because we permute for 20 times, 0.05 correspond to smaller than the smallest in the extreme value vector, if we permute 200 times 0.05 will correspond to smaller than the 10th smallest value in extreme stats vector
B = bwboundaries(binary_sig_stats);
for k = 1:length(B)
    boundary_size(k)=size(B{k},1);
end
[sort_size,I]=sort(boundary_size,'descend');
for k = 1:1
    boundary = B{I(k)};
    interp_boundx=5*(boundary(:,1)-1);
    interp_boundy=282-5*(boundary(:,2)-1);
    plot(interp_boundx, interp_boundy, 'r', 'LineWidth', 2)
end