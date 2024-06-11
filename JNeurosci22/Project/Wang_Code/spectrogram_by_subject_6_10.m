project_dir = 'F:\CCLab\Projects\LFP_project_v0.2';
output_database = 'Wang_Database';
%%
% normalizer = 3
load(fullfile(project_dir, output_database, 'groups_5_5.mat'));
load(fullfile(project_dir, output_database, 'global_params_v2.mat'));
load(fullfile(project_dir, output_database, 'cwt_post_cohen_chronux_v2.mat'));
load(fullfile(project_dir, output_database, 'cwt_pre_cohen_chronux_v2.mat'));
%%  Stage cutoff 
last_young_session = [20, 26, 66, 27];
%   Convert  dummy code
[~, last_young_session] = max(last_young_session == sessions', [], 1);
% 
target_subjects = [1,3,4];
cwt_young_ps = cell(numel(target_subjects), 2, 2);
cwt_young_as = cell(numel(target_subjects), 2, 2);
cwt_adult_ps = cell(numel(target_subjects), 2, 2);
cwt_adult_as = cell(numel(target_subjects), 2, 2);
tic
for i = 1:numel(target_subjects)
    monk_ind = target_subjects(i)
    for a = 1:size(groups, 1)
        for b = 1:size(groups, 2)
            for c = 1:size(groups, 3)
                c
    %             for d = 1:size(groups, 4)
                for d = monk_ind
                    % Loops through each kernel type
                    for kernel_switch = 1:numel(kernel_flag)
                        p_1 = squeeze(cwt_1{kernel_switch}(a, b, c, d, :, :));
                        p_2 = squeeze(cwt_2{kernel_switch}(a, b, c, d, :, :));
                        % 
                        if or(any(isnan(p_1), 'all'), any(isnan(p_2), 'all'))
                            continue
                        end
                        if c > last_young_session(d)
                            if b == 1
                                cwt_adult_ps{i, 1, kernel_switch} = cat(3, cwt_adult_ps{i, 1, kernel_switch}, p_1);
                                cwt_adult_ps{i, 2, kernel_switch} = cat(3, cwt_adult_ps{i, 2, kernel_switch}, p_2);
                            else
                                cwt_adult_as{i, 1, kernel_switch} = cat(3, cwt_adult_as{i, 1, kernel_switch}, p_1);
                                cwt_adult_as{i, 2, kernel_switch} = cat(3, cwt_adult_as{i, 2, kernel_switch}, p_2);
                            end
                        else
                            if b == 1
                                cwt_young_ps{i, 1, kernel_switch} = cat(3, cwt_young_ps{i, 1, kernel_switch}, p_1);
                                cwt_young_ps{i, 2, kernel_switch} = cat(3, cwt_young_ps{i, 2, kernel_switch}, p_2);
                            else
                                cwt_young_as{i, 1, kernel_switch} = cat(3, cwt_young_as{i, 1, kernel_switch}, p_1);
                                cwt_young_as{i, 2, kernel_switch} = cat(3, cwt_young_as{i, 2, kernel_switch}, p_2);
                            end
                        end
                    end
                end
                toc
            end
        end
    end
end
%%  Save by subject CWT
% fname_ = sprintf(...
%     'cwt_by_subject_norm_%1.0f_pre_%.1f_post_%.1f_kernel1_%d_kernel2_%d.mat', ...
%     normalizer, pre_dur, post_dur, kernel_flag(1), kernel_flag(2)...
%     );
fname_ = 'cwt_by_subject_v2.mat';
save(...
    fullfile(project_dir, output_database, fname_), ...
    'cwt_young_ps' , 'cwt_young_as' , 'cwt_adult_ps' , 'cwt_adult_as'...
    )
