function new_repo = zw_repo_cwt_temp_1_30_22(...
    repo, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    )
%%
% Temporary file used to demonstrate Fiona's data
%%  Pad to avoid edge effects
% Taper window 500 ms
pre_pad = 0.3 * fs; %  number of samples to pad towards pre-trial
post_pad = 0.3 * fs; % number pf samples length to pad towards post-trial
%%
global n_repo
n_repo = numel(repo);
new_repo = struct('class', []);
for i =1:n_repo
    new_repo(i).class = struct();
    new_repo(i).class.cue_cwt = zeros([0, 0, 0]);
    new_repo(i).class.sac_cwt = zeros([0, 0, 0]);
    new_repo(i).class.norm = zeros([0, 0]);
end
global counter;
counter = 0;

tic
counter_q = parallel.pool.DataQueue;
afterEach(counter_q, @counter_check);
% parfor i = 1:numel(repo)
for i = 1:numel(repo)
    i
    for j = 1:numel(repo(i).class)
        if isempty(repo(i).class(j).ntr)
            continue
        end
        for k  =1:numel(repo(i).class(j).ntr)
            pre_trial_T = repo(i).class(j).ntr(k).Cue_onT + 1; % Assuming 1-second fixation
            cue_sample = floor((cue_dur + 1)*fs); % Cue-aligned samples from trial onset
            pre_sample_ava =  floor(pre_trial_T * fs);
            pre_sample_req = pre_pad; % Pad backward in time
            post_sample_ava = numel(repo(i).class(j).ntr(k).LFP) - floor(pre_trial_T * fs) - 1*fs ;
            post_sample_req = post_pad + floor(cue_sample(2)); %  Pad forward in time starting from trial onset
            pre_diff = pre_sample_ava - pre_sample_req;
            if pre_diff < 0 %   Pad zeros when data is too short
                pre_sample_ava = pre_sample_req;
                repo(i).class(j).ntr(k).LFP = [zeros(-pre_diff, 1); repo(i).class(j).ntr(k).LFP];
            end
            post_diff = post_sample_ava - post_sample_req;
            if post_diff < 0
                post_sample_ava = post_sample_req;
                repo(i).class(j).ntr(k).LFP = [repo(i).class(j).ntr(k).LFP; zeros(-post_diff, 1)];
            end
                % Convert
                sig = repo(i).class(j).ntr(k).LFP';
                sig = sig - mean(sig, 2);
                % Filtering
%                 sig = zw_lfp_filt(sig);
                pre_sig = sig((pre_sample_ava - pre_sample_req + 1):(pre_sample_ava));
                post_sig = sig((pre_sample_ava + 1):(pre_sample_ava + post_sample_req));
                % Calls ZW_CWT_VER_2 for each trial
                [cue_cwt_, sac_cwt_, norm_] = zw_cwt_new(...
                    pre_sig, post_sig, normalizer, ...
                    cue_sample, cue_sample, ...
                    kernel_flag, down_sample, f_range, fs...
                    );
%             else % Trial output assigned NaN if too short
%                 cue_cwt_ = nan(numel(f_range), (diff(cue_sample) + 1)/down_sample);
%                 sac_cwt_ = nan(numel(f_range), (diff(sac_sample) + 1)/down_sample);
%                 norm_ = nan;
%             end
            new_repo(i).class(j).cue_cwt(k, :, :) = cue_cwt_;
%             new_repo(i).class(j).sac_cwt(k, :, :) = sac_cwt_;
            new_repo(i).class(j).norm(k, :) = norm_;
        end
    end
%     toc
    send(counter_q, 1);
end
end
function counter_check(q_sent)
global counter
global n_repo
counter = counter + q_sent;
if mod(counter, 100) < 1
    timer = toc
    sprintf(...
        '%d/%d, %.2f hours left', ...
        counter, ...
        n_repo, ...
        timer/(counter)*(n_repo - counter)/60/60)
end
end