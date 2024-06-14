all_target=find_index(find_sigsample_CS2);
all_target_f=f_sample(find_sigsample_CS2,:);
[temp_f,get_index]=sort(all_target_f(:,2),'descend');
final_targ_index=all_target(get_index(1:10));