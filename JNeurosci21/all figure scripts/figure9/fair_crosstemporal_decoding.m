function fair_crosstemporal_decoding(cell_group,mode)
load('allcell_post.mat');  %analyse pre or post training
cell_from_period=2; %1 from sample period, 2 from sampledelay period
find_index=allcell_info{10};
if cell_from_period==1
NMS_index=find_index(unique([allcell_info{3};allcell_info{4};allcell_info{5};allcell_info{6}]));
LS_index=find_index(unique([allcell_info{7};allcell_info{8};allcell_info{9}]));
num_select=min([length(NMS_index),length(LS_index)]);
else
NMS_index=find_index(unique([allcell_info{11};allcell_info{12};allcell_info{13};allcell_info{14}]));
LS_index=find_index(unique([allcell_info{15};allcell_info{16};allcell_info{17}]));
num_select=min([length(NMS_index),length(LS_index)]);
end
%cell_group=1;  %1 represent LS, 2 represent NMS
%mode=2; %mode1 decode cue, model decode sample, model decode match/nonmatch 

for r=1:1
if cell_group==1
temp_vec=randperm(length(LS_index));
temp_index=LS_index(temp_vec(1:num_select));
else
temp_vec=randperm(length(NMS_index));
temp_index=NMS_index(temp_vec(1:num_select));
end
remove_index=find(temp_index==1523);
temp_index(remove_index)=[];
remove_index=find(temp_index==1524);
temp_index(remove_index)=[];
remove_index=find(temp_index==255);
temp_index(remove_index)=[];
[x,y,crv_cell,population_response]=normalize_SVM_decoder(0.4,0.1,temp_index,mode);
accuracy_mat=zeros(length(x),length(x));

match_label=[ones(1,48),zeros(1,48)];
cue_label=[];
sample_label=[];
lookup_table=[5,6,7,8,1,2,3,4];

for i=1:8
cue_label=[cue_label,i*ones(1,6)];
sample_label=[sample_label,i*ones(1,6)];
end
for i=1:8
cue_label=[cue_label,i*ones(1,6)];
sample_label=[sample_label,lookup_table(i)*ones(1,6)];
end

for i=1:length(y) %loop through training time point
    disp(i);
    temp_crv=crv_cell{i};
    temp_partition=temp_crv.Partition;
    for j=1:length(y)  %loop through testing time point 
       if i==j
           accuracy_mat(i,j)=y(i);
       else
           for p=1:10  % 10 fold cross validation
               fold_decoder=temp_crv.Trained{p};
               fold_test_index=find(training(temp_partition,p)==0);
               label=predict(temp_crv.Trained{p},population_response(:,fold_test_index,j)');
               if mode==1
               correct_count(p)=length(find(label==cue_label(fold_test_index)')); %change what to decode as needed
               elseif mode==2
               correct_count(p)=length(find(label==sample_label(fold_test_index)')); %change what to decode as needed
               else
               correct_count(p)=length(find(label==match_label(fold_test_index)')); %change what to decode as needed
               end
               all_count(p)=length(fold_test_index);
           end
           accuracy_mat(i,j)=sum(correct_count)/sum(all_count);
       end
    end
end
all_results(r,:,:)=accuracy_mat;
end
%accuracy_mat=squeeze(mean(all_results));
accuracy_mat=squeeze(all_results);
savename=['post_delay_',num2str(cell_group),'_',num2str(mode),'.mat'];
%save(savename,'all_results');
[X,Y] = meshgrid(-0.8:0.1:4.8);
[Xq,Yq] = meshgrid(-0.8:0.02:4.8);
Vq = interp2(X,Y,accuracy_mat,Xq,Yq);
imagesc(flipud(Vq));
colormap parula;
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
yticklabels({'4s','3s','2s','1s','0s'});
xticks([40 90 140 190 240]);
xticklabels({'0s','1s','2s','3s','4s'});
%ylabel('training time');
%xlabel('testing time');
caxis([0.125,0.7]);
disp(size(temp_index));
end