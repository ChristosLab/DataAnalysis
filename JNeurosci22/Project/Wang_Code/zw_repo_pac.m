function new_repo = zw_repo_pac(repo, cfg, target_lfp)
fs = ft_getopt(cfg,'Fs',[]);
durs = ft_getopt(cfg,'durs',[]);
phase_freqs = ft_getopt(cfg,'phase_freqs',[]);
amp_freqs = ft_getopt(cfg,'amp_freqs',[max(phase_freqs):2:fs/2]);
%%
global n_target
n_repo = numel(repo);
n_target = numel(target_lfp);
new_repo = struct('class', []);
for i =1:n_repo
    new_repo(i).class = struct();
    new_repo(i).class.pac = zeros([0, 0, 0]);
    new_repo(i).class.p = zeros([0, 0, 0]);
    new_repo(i).class.miu = zeros([0, 0, 0]);
    new_repo(i).class.delta = zeros([0, 0, 0]);
end
global counter;
counter = 0;

tic
counter_q = parallel.pool.DataQueue;
afterEach(counter_q, @counter_check);
parfor i = 1:numel(repo)
    if ~any(target_lfp == i)
        continue
    end
    i
    for j = 1:numel(repo(i).class)
        if isempty(repo(i).class(j).ntr)
            continue
        end
%       
        sig_window = [-1 + 1/fs, 3];
        sig = zeros(0, diff(sig_window) * fs + 1);
        for k  =1:numel(repo(i).class(j).ntr)
            if isempty(repo(i).class(j).ntr(k).LFP)
                continue
            end
            Cue_onT = repo(i).class(j).ntr(k).Cue_onT;
            sig_sample = floor((sig_window(1) + Cue_onT)* fs) + [0, size(sig, 2) - 1];
            sig = [sig; repo(i).class(j).ntr(k).LFP(sig_sample(1):sig_sample(2))'];
        end
                cfg_ = cfg;
                cfg_.mask = floor((durs - sig_window(1)) * fs);
                [pac_, miu_, delta_] = PACmeg_surrogate(cfg_, sig);
%                 pac_ = zeros([0, 0, 0, 0, 0]);
            new_repo(i).class(j).pac(:, :, :) = pac_;
            new_repo(i).class(j).miu(:, :, :) = miu_;
            new_repo(i).class(j).delta(:, :, :) = delta_;
%         for k  =1:numel(repo(i).class(j).ntr)
%             if isempty(repo(i).class(j).ntr(k).LFP)
%                 pac_ = nan(size(mask, 1), numel(amp_freqs), numel(phase_freqs));
%                 new_repo(i).class(j).pac_(k, :, :, :) = pac_;
%                 continue
%             end
%             Cue_onT = repo(i).class(j).ntr(k).Cue_onT;
%                 sig = repo(i).class(j).ntr(k).LFP';
%                 % Calls PAC_meg for each trial
%                 cfg_ = cfg;
%                 cfg_.mask = floor((durs + Cue_onT) * fs);
%                 pac_ = PACmeg(cfg_, sig);
%             new_repo(i).class(j).pac(k, :, :, :) = pac_;
%         end
    end
%     toc
    send(counter_q, 1);
end
end
function counter_check(q_sent)
global counter
global n_target
counter = counter + q_sent;
if mod(counter, 5) < 1
    timer = toc
    sprintf(...
        '%d/%d, %.2f hours left', ...
        counter, ...
        n_target, ...
        timer/(counter)*(n_target - counter)/60/60)
end
end