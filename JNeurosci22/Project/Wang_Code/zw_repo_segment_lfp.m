function new_repo = zw_repo_segment_lfp(...
    repo, ...
    cue_dur, fs...
    )
% Taper window 500 ms
%%
load('D:\Database\Zhengyang_Wang\Projects\LFP_project_v0.2\Wang_Database\lfp_tbl.mat')
target_lfp = find(lfp_tbl.task_id == 1);
global n_repo
n_repo = numel(repo);
new_repo = struct('class', []);
for i =1:n_repo
    new_repo(i).class = struct();
    new_repo(i).class.lfp_cue = zeros([0, 0, 0]);
end
global counter;
counter = 0;

tic
counter_q = parallel.pool.DataQueue;
afterEach(counter_q, @counter_check);
for i = 1:numel(repo)
%     i
    if ~any(target_lfp == i)
        continue
    end
    for j = 1:numel(repo(i).class)
        if isempty(repo(i).class(j).ntr)
            continue
        end
%         if ~isfield(repo(i).class(j).ntr, 'Saccade_onT')
%             continue
%         end
        cue_sig = nan(numel(repo(i).class(j).ntr), round(diff(cue_dur)*fs));
        for k  =1:numel(repo(i).class(j).ntr)
            if isempty(repo(i).class(j).ntr(k).LFP)
                continue
            end
            %             if isempty(repo(i).class(j).ntr(k).Saccade_onT) || isnan(repo(i).class(j).ntr(k).Saccade_onT)
            %                 continue
            %             end
            cue_sample = floor(cue_dur*fs) + floor(repo(i).class(j).ntr(k).Cue_onT*fs);
            %             sac_sample = floor((sac_dur + repo(i).class(j).ntr(k).Saccade_onT)*fs);
            %
            %
            % Convert
            cue_sig(k, :) = repo(i).class(j).ntr(k).LFP((cue_sample(1) + 1):cue_sample(end))';
        end
        % Filtering
        cue_sig = zw_lfp_filt(cue_sig);
        new_repo(i).class(j).lfp_cue = cue_sig;
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