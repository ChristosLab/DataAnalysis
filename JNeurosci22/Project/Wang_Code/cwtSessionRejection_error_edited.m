zw_setpath
fname_ = 'lfp_tbl.mat';
load(fullfile(project_dir, output_database, fname_), 'lfp_tbl');
%%
fname_ = 'cwt_error_repo.mat';
load(fullfile(project_dir, output_database, fname_), 'cwt_error_repo');
t_y = [t{1, :}];
t_a = [t{2, :}];
t_ya = {t_y, t_a};
error_session_num_old = [];
for i = 1:numel(t_ya)
    for j = 1:numel(t_ya{i})
        current_session = t_ya{i}(j);
        if isempty(cwt_error_repo(current_session).class(1).norm) && numel(cwt_error_repo(current_session).class) == 1
            continue
        end
        error_session_num_old = [error_session_num_old, current_session];
    end
end
% 
fname_ = 'cwt_error_repo_edited.mat';
load(fullfile(project_dir, output_database, fname_), 'cwt_error_repo');
error_session_num_new = [];
for i = 1:numel(t_ya)
    for j = 1:numel(t_ya{i})
        current_session = t_ya{i}(j);
        if isempty(cwt_error_repo(current_session).class(1).norm) && numel(cwt_error_repo(current_session).class) == 1
            continue
        end
        error_session_num_new = [error_session_num_new, current_session];
    end
end
%%
global list_i flag_list ind_list
%%
target_cwts = setxor(error_session_num_new, error_session_num_old);
%%
ind_list = zeros(0, 2);
for i =  target_cwts
    for j = 1:numel(cwt_error_repo(i).class)
        if ~isempty(cwt_error_repo(i).class(j).cue_cwt)
            ind_list = [ind_list(); [i, j]];
        end
    end
    i
end
%%
global log_dir
log_dir = fullfile(project_dir, output_database, 'error_artifact_rejection_list_additional.mat');
%%  Initiate values
if isfile(log_dir)
    load(log_dir)
else
    flag_list = zeros(size(ind_list, 1), 1);
    list_i = 1;
end
%%  Run GUI to reject trials and store flags in  flag_list
f = figure('Units', 'pixels', 'Position', [25, 25, 900, 900]);
% ax = axes(f);s
% ax.Units = 'pixels';
% ax.Position = [75 75 600 300];
c_flag = uicontrol('Style', 'radiobutton', 'Position', [600, 25, 100, 50]);
c_flag.String = {'Rejected!'};
c_flag.Callback = @radioToggled;
c_forward = uicontrol('Position', [425, 25, 100, 50]);
c_forward.String = 'Next trial';
c_forward.Callback = {@forwardButtonPushed, cwt_error_repo, c_flag};
c_backward = uicontrol('Position', [225, 25, 100, 50]);
c_backward.String = 'Previous trial';
c_backward.Callback = {@backwardButtonPushed, cwt_error_repo, c_flag};
f.KeyPressFcn = {@keyPress, c_backward.Callback, c_forward.Callback, c_flag};
%   Inititate first plot
cwt = cwt_error_repo(ind_list(list_i, 1)).class(ind_list(list_i, 2)).cue_cwt;
plot_cwt(cwt)
%%  Reject sessions
k = 0;
session_flag_list = zeros(size(cwt_error_repo));
for i = 1:numel(target_cwts)
    for j = 1:numel(cwt_error_repo(target_cwts(i)).class)
        if isempty(cwt_error_repo(target_cwts(i)).class(j).norm)
            continue
        end
            k = k + 1;
            if flag_list(k)
                session_flag_list(target_cwts(i)) = 1;
            end
    end
    i
end
%%  Printing current UI
fig_name_ = sprintf('%s Class %d', lfp_tbl.Filename{ind_list(list_i, 1)}, ind_list(list_i, 2));
sgtitle(fig_name_, 'FontWeight', 'bold', 'Interpreter', 'none')
print(gcf, fullfile(project_dir, fig_lib, fig_name_), '-dpdf', '-fillpage');
%%
fname_ = 'maunal_session_rejection_flag_error_additional.mat';
save(fullfile(project_dir, output_database, fname_), 'session_flag_list', 'flag_list', 'ind_list', 'target_cwts');
%%
fname_ = 'maunal_session_rejection_flag_error.mat';
load(fullfile(project_dir, output_database, fname_), 'session_flag_list', 'flag_list', 'ind_list', 'target_cwts');
%%  GUI functions
function keyPress(src, event, varargin)
switch event.Key
    case 'leftarrow'
        varargin{1}{1}(src, [], varargin{1}{2:end});
    case 'rightarrow'
        varargin{2}{1}(src, [], varargin{2}{2:end});
    case 'space'
        varargin{3}.Value = ~varargin{3}.Value;
end
end

function forwardButtonPushed(src,event,cwt_error_repo, c_flag)
global list_i flag_list ind_list
if list_i == size(ind_list, 1)
    disp('Finished')
    return
end
flag_list(list_i) = c_flag.Value;
list_i = list_i + 1;
cwt = cwt_error_repo(ind_list(list_i, 1)).class(ind_list(list_i, 2)).cue_cwt;
delete(findall(src, 'type', 'axes'));
delete(findall(src, 'type', 'text'));
try 
plot_cwt(cwt)
c_flag.Value = flag_list(list_i);
log_flag_list();
catch
    forwardButtonPushed(src,event,cwt_error_repo, c_flag)
end
end

function backwardButtonPushed(src,event,cwt_error_repo, c_flag)
global list_i flag_list ind_list
if list_i == 1
    disp('No more previous trials')
    return
end
flag_list(list_i) = c_flag.Value;
list_i = list_i - 1;
cwt = cwt_error_repo(ind_list(list_i, 1)).class(ind_list(list_i, 2)).cue_cwt;
delete(findall(src, 'type', 'axes'));
delete(findall(src, 'type', 'text'));
plot_cwt(cwt)
c_flag.Value = flag_list(list_i);
log_flag_list();
end
% 
% function radioUpdated(src, event)
% global list_i flag_list ind_list
% end
function radioToggled(hObject, event)
global list_i flag_list
flag_list(list_i) = hObject.Value;
log_flag_list();
end
function log_flag_list()
global flag_list list_i log_dir
save(log_dir, 'flag_list', 'list_i');
end
function plot_cwt(cwt)
global list_i ind_list
set(gcf, 'Name',sprintf('Session %d Class %d', ind_list(list_i, 1), ind_list(list_i, 2)))
n_plot_per_row = 4;
n_trial = size(cwt, 1);
for i = 1:n_trial
    subplot(2*ceil(n_trial/n_plot_per_row), n_plot_per_row, mod(i -1, n_plot_per_row) + 1 + 2*n_plot_per_row*floor((i-1)/n_plot_per_row))
    imagesc([-1, 3], [2, 128], squeeze(cwt(i, :, :))./squeeze(mean(cwt(i, :,26:50), 3))');
    set(gca, 'YDir', 'normal')
    subplot(2*ceil(n_trial/n_plot_per_row), n_plot_per_row, mod(i -1, n_plot_per_row) + 1 + 2*n_plot_per_row*floor((i-1)/n_plot_per_row) + n_plot_per_row)
    imagesc([-1, 3], [2, 128], squeeze(cwt(i, :, :)));
    set(gca, 'YDir', 'normal')
end

% sgtitle(sprintf('Session %d Class %d', ind_list(list_i, 1), ind_list(list_i, 1)));
end
%%