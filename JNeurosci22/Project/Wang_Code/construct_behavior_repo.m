%Construct behavior repo matching the LFP repo
behavior_repo = struct([]);
y_out_counter = 0;
for i = [site_mod_cat{:}, site_mod_cat_nonsig{:}]
    load(fullfile(mat_database, [lfp_tbl.Filename{i}, '.mat']));
    behavior_repo(i).version    = AllData.version;
    behavior_repo(i).parameters = AllData.parameters;
    behavior_repo(i).Filename   = lfp_tbl.Filename{i};
    [n_trial, s_rate, y_out] = compute_performance(AllData.trials);
    y_out_counter = y_out_counter + y_out;
    behavior_repo(i).n_trial    = n_trial;
    behavior_repo(i).s_rate     = s_rate;
    clear AllData
end
%%
fname_ = 'odr_behavior_repo';
save(fullfile(project_dir, output_database, fname_), 'behavior_repo');
%%
for i = 1:numel()
    load(fullfile(mat_database, [lfp_tbl.Filename{i}, '.mat']));
    behavior_repo(i).version    = AllData.version;
    behavior_repo(i).parameters = AllData.parameters;
    behavior_repo(i).Filename   = lfp_tbl.Filename{i};
    [n_trial, s_rate] = compute_performance(AllData.trials);
    behavior_repo(i).n_trial    = n_trial;
    behavior_repo(i).s_rate     = s_rate;
    clear AllData
end

%% gamma delay - behavior
figure
plot(mean([lfp_mod_cue([site_mod_cat{1, :, :}]', 76:150, 3)], 2), [behavior_repo([site_mod_cat{1, :, :}]).s_rate], 'r*')
hold on
plot(mean([lfp_mod_cue([site_mod_cat{2, :, :}]', 76:150, 3)], 2), [behavior_repo([site_mod_cat{2, :, :}]).s_rate], 'b*')
%% gamma delay - behavior (tuned sites only only)
figure
plot(mean([lfp_mod_cue([site_tuning_cat{1, :, :, 2}]', 76:150, 3)], 2), [behavior_repo([site_tuning_cat{1, :, :, 2}]).s_rate], 'r*')
hold on
plot(mean([lfp_mod_cue([site_tuning_cat{2, :, :, 2}]', 76:150, 3)], 2), [behavior_repo([site_tuning_cat{2, :, :, 2}]).s_rate], 'b*')
%% high gamma delay - behavior
figure
plot(mean([lfp_mod_cue([site_mod_cat{1, :, :}]', 76:150, 4)], 2), [behavior_repo([site_mod_cat{1, :, :}]).s_rate], 'r*')
hold on
plot(mean([lfp_mod_cue([site_mod_cat{2, :, :}]', 76:150, 4)], 2), [behavior_repo([site_mod_cat{2, :, :}]).s_rate], 'b*')
%% gamma cue - behavior
figure
plot(mean([lfp_mod_cue([site_mod_cat{1, :, 1:2}]', 51:75, 3)], 2), [behavior_repo([site_mod_cat{1, :, 1:2}]).s_rate], 'r*')
hold on
plot(mean([lfp_mod_cue([site_mod_cat{2, :, 1:2}]', 51:75, 3)], 2), [behavior_repo([site_mod_cat{2, :, 1:2}]).s_rate], 'b*')
xl = xlim;
plot(xl, [0, 0] + mean([behavior_repo([site_mod_cat{1, :, 1:2}]).s_rate]), 'r')
plot(xl, [0, 0] + mean([behavior_repo([site_mod_cat{2, :, 1:2}]).s_rate]), 'b')
%% high gamma cue - behavior
figure
plot(mean([lfp_mod_cue([site_mod_cat{1, :, 1:2}]', 51:75, 4)], 2), [behavior_repo([site_mod_cat{1, :, 1:2}]).s_rate], 'r*')
hold on
plot(mean([lfp_mod_cue([site_mod_cat{2, :, 1:2}]', 51:75, 4)], 2), [behavior_repo([site_mod_cat{2, :, 1:2}]).s_rate], 'b*')

%%
function [n_trial, s_rate, y_out] = compute_performance(trials)
trial_counter   = 0;
s_trial_counter = 0;
temp_counter = 0;
for i = 1:numel(trials)
%     if ~isempty(trials(i).TargetIn)
%         trial_counter = trial_counter + 1;
%     end
    if ismember(trials(i).Statecode, [5,6])
        trial_counter = trial_counter + 1;
    end
    if trials(i).Statecode == 5
        temp_counter = temp_counter  + 1;
    end
    if strcmp(trials(i).Reward, 'Yes')
        s_trial_counter = s_trial_counter + 1;
    end
    n_trial = trial_counter;
    s_rate = s_trial_counter/n_trial;
end
if s_trial_counter + temp_counter == n_trial
    y_out = 1;
else 
    y_out = 0;
end
end
