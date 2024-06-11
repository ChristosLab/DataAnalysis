% repo_i = 3;
% class_i = 7;
% for i = 1:size(cwt_repo(repo_i).class(class_i).cue_cwt, 1)
%     imagesc(squeeze(cwt_repo(repo_i).class(class_i).cue_cwt(i,:,:))./squeeze(mean(cwt_repo(repo_i).class(class_i).cue_cwt(i,:,26:50), 3))');
%     set(gca, 'YDir', 'normal')
%     pause;
%     close all
% end
% %%
% plot(1,1)
% o_1 = questdlg('Is this trial an artifact', 'Artifact check', 'yes', 'no', 'no');
% if isempty(o_1)
%     o_1 = 0;
% elseif strcmp(o_1, 'yes')
%     o_1 = 1;
% elseif strcmp(o_1, 'no')
%     o_1 = 0;
% end
% o_1
%%
global list_i flag_list ind_list

%%
ind_list = zeros(0, 3);
for i = 1:numel(cwt_repo)
    for j = 1:numel(cwt_repo(i).class)
        for k = 1:size(cwt_repo(i).class(j).cue_cwt, 1)
            ind_list = [ind_list(); [i, j, k]];
        end
    end
    i
end
%%
flag_list = zeros(size(ind_list, 1), 1);
list_i = 1;
f = figure('Units', 'pixels', 'Position', [25, 25, 800, 800]);
ax = axes(f);
ax.Units = 'pixels';
ax.Position = [75 75 600 300];
% update_callback = {@, };
c_flag = uicontrol('Style', 'radiobutton', 'Position', [200, 25, 50, 50]);
c_flag.String = {'Rejected!'};
c_flag.Callback = @radioToggled;
c_forward = uicontrol('Position', [125, 25, 50, 50]);
c_forward.String = 'Next trial';
c_forward.Callback = {@forwardButtonPushed, cwt_repo, c_flag};
c_backward = uicontrol('Position', [25, 25, 50, 50]);
c_backward.String = 'Previous trial';
c_backward.Callback = {@backwardButtonPushed, cwt_repo, c_flag};
f.KeyPressFcn = {@keyPress, c_backward.Callback, c_forward.Callback, c_flag};

cwt = squeeze(cwt_repo(ind_list(list_i, 1)).class(ind_list(list_i, 2)).cue_cwt(ind_list(list_i, 3),:,:));
plot_cwt(cwt)

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

function forwardButtonPushed(src,event,cwt_repo, c_flag)
global list_i flag_list ind_list
if list_i == size(ind_list, 1)
    disp('Finished')
    return
end
flag_list(list_i) = c_flag.Value;
list_i = list_i + 1;
ind = ind_list(list_i, :);
cwt = squeeze(cwt_repo(ind_list(list_i, 1)).class(ind_list(list_i, 2)).cue_cwt(ind_list(list_i, 3),:,:));
plot_cwt(cwt)
c_flag.Value = flag_list(list_i);
end

function backwardButtonPushed(src,event,cwt_repo, c_flag)
global list_i flag_list ind_list
if list_i == 1
    disp('No more previous trials')
    return
end
flag_list(list_i) = c_flag.Value;
list_i = list_i - 1;
cwt = squeeze(cwt_repo(ind_list(list_i, 1)).class(ind_list(list_i, 2)).cue_cwt(ind_list(list_i, 3),:,:));
plot_cwt(cwt)
c_flag.Value = flag_list(list_i);
end
% 
% function radioUpdated(src, event)
% global list_i flag_list ind_list
% end
function radioToggled(hObject, event)
global list_i flag_list
flag_list(list_i) = hObject.Value;
end
function plot_cwt(cwt)
subplot(2,1,1)
imagesc(cwt./mean(cwt(:,26:50), 2));
set(gca, 'YDir', 'normal')
subplot(2,1,2)
imagesc(cwt);
set(gca, 'YDir', 'normal')
end