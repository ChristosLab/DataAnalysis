%%  Stage cutoff 
last_young_session = [20, 26, 66, 27];
%   Convert  dummy code
[~, last_young_session] = max(last_young_session == sessions', [], 1);
%%  Load stft
interval_pre = [-1.5, 0];
interval_post = [-1, 2.5];
interval_pre_for_pow = [-1, 0];
window = 250;
overlap = 240;
f_range = 0:2:100;
fname_ = sprintf('stft_pre_norm_%.1f_%.1f_%.1f_%.1f.mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2));
load(fullfile(project_dir, output_database, fname_), 'stft_1', 'window', 'overlap', 'f_range', 'interval_pre', 'interval_post');

fname_ = sprintf('stft_tot_norm_%.1f_%.1f_%.1f_%.1f.mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2));
load(fullfile(project_dir, output_database, fname_), 'stft_2', 'window', 'overlap', 'f_range', 'interval_pre', 'interval_post');

fname_ = sprintf('rms_pow_%.1f_%.1f_%.1f_%.1f.mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2));
load(fullfile(project_dir, output_database, fname_), 'abs_pow', 'interval_pre', 'interval_post');

%%
normalizer = 1;
%%  Aggragate STFT
if normalizer == 0
    stft = stft_1;
elseif normalizer == 1
    stft = stft_2;
end

target_subjects = [1,3,4];
stft_young_ps = cell(numel(target_subjects), 1);
stft_young_as = cell(numel(target_subjects), 1);
stft_adult_ps = cell(numel(target_subjects), 1);
stft_adult_as = cell(numel(target_subjects), 1);
tic
for i = 1:numel(target_subjects)
    monk_ind = target_subjects(i)
    for a = 1:size(groups, 1)
        for b = 1:size(groups, 2)
            for c = 1:size(groups, 3)
                c
    %             for d = 1:size(groups, 4)
                for d = monk_ind
                    p_ = squeeze(stft(a, b, c, d, :, :));
                    if any(isnan(p_))
                        continue
                    end
                    if c > last_young_session(d)
                        if b == 1
                            stft_adult_ps{i} = cat(3, stft_adult_ps{i}, p_);
                        else
                            stft_adult_as{i} = cat(3, stft_adult_as{i}, p_);
                        end
                    else
                        if b == 1
                            stft_young_ps{i} = cat(3, stft_young_ps{i}, p_);
                        else
                            stft_young_as{i} = cat(3, stft_young_as{i}, p_);
                        end
                    end
                end
                toc
            end
        end
    end
end
%%  Save by subject STFT
fname_ = sprintf(...
    'stft_by_subject_norm_%1.0f_%.1f_%.1f_%.1f_%.1f.mat', ...
    normalizer, interval_pre(1),interval_pre(2),interval_post(1),interval_post(2)...
    );
save(...
    fullfile(project_dir, output_database, fname_), ...
    'stft_young_ps' , 'stft_young_as' , 'stft_adult_ps' , 'stft_adult_as'...
    )
%%
%%  Aggregate CWT
target_subjects = [1,3,4];
cwt_young_ps = cell(numel(target_subjects), 2);
cwt_young_as = cell(numel(target_subjects), 2);
cwt_adult_ps = cell(numel(target_subjects), 2);
cwt_adult_as = cell(numel(target_subjects), 2);
tic
for i = 1:numel(target_subjects)
    monk_ind = target_subjects(i)
    for a = 1:size(groups, 1)
        for b = 1:size(groups, 2)
            for c = 1:size(groups, 3)
                c
    %             for d = 1:size(groups, 4)
                for d = monk_ind
                    p_1 = squeeze(cwt_1(a, b, c, d, :, :));
                    p_2 = squeeze(cwt_2(a, b, c, d, :, :));
                    if or(any(isnan(p_1), 'all'), any(isnan(p_2), 'all'))
                        continue
                    end
                    if c > last_young_session(d)
                        if b == 1
                            cwt_adult_ps{i, 1} = cat(3, cwt_adult_ps{i, 1}, p_1);
                            cwt_adult_ps{i, 2} = cat(3, cwt_adult_ps{i, 2}, p_2);
                        else
                            cwt_adult_as{i, 1} = cat(3, cwt_adult_as{i, 1}, p_1);
                            cwt_adult_as{i, 2} = cat(3, cwt_adult_as{i, 2}, p_2);
                        end
                    else
                        if b == 1
                            cwt_young_ps{i, 1} = cat(3, cwt_young_ps{i, 1}, p_1);
                            cwt_young_ps{i, 2} = cat(3, cwt_young_ps{i, 2}, p_2);
                        else
                            cwt_young_as{i, 1} = cat(3, cwt_young_as{i, 1}, p_1);
                            cwt_young_as{i, 2} = cat(3, cwt_young_as{i, 2}, p_2);
                        end
                    end
                end
                toc
            end
        end
    end
end
%%  Save by subject CWT
normalizer = 0;
fname_ = sprintf(...
    'cwt_by_subject_norm_%1.0f_%.1f_%.1f_%.1f_%.1f.mat', ...
    normalizer, interval_pre(1),interval_pre(2),interval_post(1),interval_post(2)...
    );
save(...
    fullfile(project_dir, output_database, fname_), ...
    'cwt_young_ps' , 'cwt_young_as' , 'cwt_adult_ps' , 'cwt_adult_as'...
    )