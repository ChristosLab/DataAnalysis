num_behavior_thred=4;
behavior_trials_perblock=40;
match_list=[15,13,11,9,17,1,3,5,7]; %match class
nonmatch_list=[16,14,12,10,2,4,6,8];%nonmatch class
load('all_behavior_info.mat');
behavior_session=behavior_info{1};
behavior_type=behavior_info{2};
behavior_statecode=behavior_info{4};
behavior_reward=behavior_info{5};
behavior_class=behavior_info{6};
find_msng=find(behavior_type==2);

for i=1:length(find_msng)

    cell_statecode=behavior_statecode(find_msng(i));
    cell_statecode=cell2mat(cell_statecode{1});
    cell_trialclass=cell2mat(behavior_class{find_msng(i)});
    cell_reward=behavior_reward(find_msng(i));
    cell_reward=cell_reward{1};
    cell_reward=strcmp(cell_reward,'Yes');
    finished_index=find(cell_statecode>5); %index for finished trials in all trials
    correct_infinish=cell_reward(find(cell_statecode>5));%correct trials in all finished trials
    cell_trialclass=cell_trialclass(find(cell_statecode>5));%class in all finished trials

    cell_behavior_blocks(i)=floor(length(finished_index)/behavior_trials_perblock);
    cumulative_correct=cumsum(correct_infinish);
    for c=1:cell_behavior_blocks(i)
        correct_inblock=sum(correct_infinish(behavior_trials_perblock*(c-1)+1:behavior_trials_perblock*c));
        correct_proportion(i,c)=correct_inblock/behavior_trials_perblock;
        for j=1:17
            finished_class_trials=[correct_infinish & (cell_trialclass==j)];
            class_correct_inbolock=sum(finished_class_trials(behavior_trials_perblock*(c-1)+1:behavior_trials_perblock*c));
            correct_class_proportion(i,c,j)=class_correct_inbolock/sum(cell_trialclass(behavior_trials_perblock*(c-1)+1:behavior_trials_perblock*c)==j);
        end
    end  
    
end
qualify_session=find(cell_behavior_blocks>=num_behavior_thred);
figure;
plot(nanmean(correct_proportion(qualify_session,1:num_behavior_thred)));
title('mean peform, 160 finished trials');

figure;
mean_class_proportion=squeeze(nanmean(correct_class_proportion(qualify_session,1:num_behavior_thred,:)));
for i=1:num_behavior_thred
nonmatch_mat(i,:)=[mean(mean_class_proportion(i,nonmatch_list([4,5]))),mean(mean_class_proportion(i,nonmatch_list([3,6]))),mean(mean_class_proportion(i,nonmatch_list([2,7]))),mean(mean_class_proportion(i,nonmatch_list([1,8])))];
match_mat(i,:)=[mean(mean_class_proportion(i,match_list([4,6]))),mean(mean_class_proportion(i,match_list([3,7]))),mean(mean_class_proportion(i,match_list([2,8]))),mean(mean_class_proportion(i,match_list([1,9]))),mean_class_proportion(i,match_list(5))];
end
b=bar(nonmatch_mat*100);
b(1).FaceColor=[230,230,230]/255;
b(2).FaceColor=[179,179,179]/255;
b(3).FaceColor=[128,128,128]/255;
b(4).FaceColor=[77,77,77]/255;
%title("Behavior performance,nonmatch,40 finished trial non overlap block" );
%xlabel('block order');
%ylabel('Percent correct');
ylim([0,140]);
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2);
yticks([0,20,40,60,80,100]);
box off;
lgd=legend('11 deg','22 deg','45 deg','90 deg');
lgd.FontSize = 14;
legend boxoff;

figure;
b=bar(match_mat*100);
b(1).FaceColor=[230,230,230]/255;
b(2).FaceColor=[179,179,179]/255;
b(3).FaceColor=[128,128,128]/255;
b(4).FaceColor=[77,77,77]/255;
b(5).FaceColor=[26,26,26]/255;
%title("Behavior performance,match,40 finished trial non overlap block" );
%xlabel('block order');
%ylabel('Percent correct');
lgd=legend('11 deg','22 deg','45 deg','90 deg','0 deg');
ylim([0,140]);
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2);
yticks([0,20,40,60,80,100]);
box off;
lgd=legend('11 deg','22 deg','45 deg','90 deg','0 deg');
lgd.FontSize = 14;
legend boxoff;
