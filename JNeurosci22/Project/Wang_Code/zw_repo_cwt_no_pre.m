function new_repo = zw_repo_cwt_no_pre(...
    repo, ...
    post_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    )
%%  Pad to avoid edge effects
pre_pad = 0.1 * fs; %  number of samples to pad towards pre-trial
post_pad = 1 * fs; % number pf samples length to pad towards post-trial
%%
global n_repo
n_repo = numel(repo);
new_repo = struct('cwt', []);
new_repo(n_repo).cwt = [];
global counter;
counter = 0;

tic
counter_q = parallel.pool.DataQueue;
afterEach(counter_q, @counter_check);
parfor i = 1:numel(repo)
    % Check if both the current and the previous trial are long enough
    pre_trial_T = repo(i).Cue_onT - 1;
    pre_sample_ava =  floor(pre_trial_T * fs);
    pre_sample_req = pre_pad; % Pad backward in time
    post_sample_ava = numel(repo(i).LFP) - floor(pre_trial_T * fs) ;
    post_sample_req = post_pad + post_dur*fs; %  Pad forward in time
    if and(...
            pre_sample_ava > pre_sample_req, ...
            post_sample_ava > post_sample_req...
            )
        % Convert 
        sig = repo(i).LFP';
        % Filtering
        sig = zw_lfp_filt(sig);
        pre_sig = sig((pre_sample_ava - pre_sample_req + 1):(pre_sample_ava));
        post_sig = sig((pre_sample_ava + 1):(pre_sample_ava + post_sample_req));
        % Calls ZW_CWT_VER_2 for each trial
        post_cwt_ = zw_cwt_no_pre(...
            pre_sig, post_sig, normalizer, ...
            pre_pad, post_pad, ...
            kernel_flag, down_sample, f_range, fs...
            );
    else % Trial output assigned NaN if too short
        post_cwt_ = nan;
    end
    new_repo(i).cwt = post_cwt_;
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