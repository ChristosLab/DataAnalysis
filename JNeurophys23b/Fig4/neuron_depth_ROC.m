% Time-resolved ROC analysis of PFC neurons at diff layers in ODR task
% J Zhu, 20230615
%% load data
clearvars
load('sig_odr_data_depth_20230906_raw_max.mat');
file_name = string(neuron_info.Filename);
neuron_info.ID = extract(file_name,1);
%% seg data/label groups
neuron_info.group(neuron_info.Depth<=800)=1;
neuron_info.group(neuron_info.Depth>800&neuron_info.Depth<=1200)=2;
neuron_info.group(neuron_info.Depth>1200)=3;
odr_data = odr_data(~isnan(neuron_info.Depth),:);
neuron_info = neuron_info(~isnan(neuron_info.Depth),:);
%% select neurons
odr_data = odr_data(neuron_info.sulcus==0,:);
neuron_info = neuron_info(neuron_info.sulcus==0,:);
%% oppo locations
opp_index = [5 6 7 8 1 2 3 4];
neuron_info.opp_del= opp_index(neuron_info.del)';
neuron_info.opp_cue = opp_index(neuron_info.cue)';
neuron_info.opp_sac= opp_index(neuron_info.sac)';
%% ROC on best cue location/opp location
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
wid = 50; % ms
noedge = 1;
for g = 1:size(group,1)
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_class = neuron_info.cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Opp_class = neuron_info.opp_cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        neuron_ROC = [];
        for n = 1:length(Neurons)
            try
                rate_cue = [];
                rate_opp_cue = [];
                [spiketrain_temp1, ntrs_temp1] = Get_spiketrain_partial_aligncue(odr_data_group(n,:),Best_class(n),[]);
                [spiketrain_temp2, ntrs_temp2] = Get_spiketrain_partial_aligncue(odr_data_group(n,:),Opp_class(n),[]);
                for ntr = 1:min(ntrs_temp1,ntrs_temp2)
                    [spks1, tlo1, thi1] = spkmtx(spiketrain_temp1(ntr),0,[-1000,3500]);
                    t1 = tlo1:thi1;
                    rate1 = 1000*smooth_es(spks1', wid, noedge)';
                    rate_cue(ntr,:) = rate1;
                    [spks2, tlo2, thi2] = spkmtx(spiketrain_temp2(ntr),0,[-1000,3500]);
                    t2 = tlo2:thi2;
                    rate2 = 1000*smooth_es(spks2', wid, noedge)';
                    rate_opp_cue(ntr,:) = rate2;
                end
                for nt = 1:size(rate_cue,2)
                    neuron_ROC(n,nt) = arrayROC(rate_cue(:,nt),rate_opp_cue(:,nt));
                end
            catch
                disp(['error processing neuron  ', num2str(Neurons(n)) '  Dir1=' num2str(Best_class(n))])
            end
        end
        NN1 = size(neuron_ROC,1);
        mr1 = mean(neuron_ROC, 1); % mean rate
        semr1 = std(neuron_ROC, 1)/sqrt(NN1); % SEM
        plt_save.t{g,s} = t1;
        plt_save.Nt{g,s} = NN1;
        plt_save.mr{g,s} = mr1;
        plt_save.semr{g,s} = semr1;
    end
end
%% plot time resolved auROC
% define a few of colors
my_color = linspecer(3);
my_color = my_color([1:3],:);
fig = figure;
clf
set(gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.85,0.35] );
for g = 1:size(group,1)
    subplot(1,3,g)
    for s = 1:size(stage,1)
        xlo = min(plt_save.t{g,s});
        xhi = max(plt_save.t{g,s});
        ylo = 0;
        yhi = max(plt_save.mr{g,s})*1.1;
        hold on;
        set(gca, 'tickdir', 'out')

        % plot errors as a shaded area
        y11 = plt_save.mr{g,s} + plt_save.semr{g,s};
        y12 = plt_save.mr{g,s} - plt_save.semr{g,s};
        fill([plt_save.t{g,s} plt_save.t{g,s}(end:-1:1)], [y11 y12(end:-1:1)], 'w', ...
            'EdgeColor', my_color(g,:)/s, 'FaceColor', my_color(g,:)/s,'FaceAlpha',0.3,'EdgeAlpha',0.3)

        % plot mean on top
        plot(plt_save.t{g,s}, plt_save.mr{g,s}, 'Color',my_color(g,:)/s, 'linewidth', 1)
        plot([0 0], [0.3 0.9], "Color",my_color(g,:)/s)
        plot([2000 2000], [0.3 0.9], "Color",my_color(g,:)/s)
    end
    axis([-1000 3500 0.3 0.9])
    xlabel('Time (ms)')
    ylabel('auROC')
    title('auROC \pm 1 SE')
    legend({'adult','','','','young'})
    %     annotation('textbox','String', "evoke best cue")
end
origUnits = fig.Units;
fig.Units = fig.PaperUnits;
fig.PaperSize = fig.Position(3:4);% set the Page Size (figure's PaperSize) to match the figure size in the Paper units
fig.Units = origUnits;% restore the original figure Units
%% plot percentage of neurons reaching different levels of ROC heatmap
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
wid = 200; % ms
noedge = 1;
figure('Position', [100 100 450 800])
for g = 3
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_class = neuron_info.cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Opp_class = neuron_info.opp_cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        neuron_ROC = [];
        for n = 1:length(Neurons)
            try
                rate_cue = [];
                rate_opp_cue = [];
                [spiketrain_temp1, ntrs_temp1] = Get_spiketrain_partial_aligncue(odr_data_group(n,:),Best_class(n),[]);
                [spiketrain_temp2, ntrs_temp2] = Get_spiketrain_partial_aligncue(odr_data_group(n,:),Opp_class(n),[]);
                for ntr = 1:min(ntrs_temp1,ntrs_temp2)
                    [spks1, tlo1, thi1] = spkmtx(spiketrain_temp1(ntr),0,[-1000,3500]);
                    t1 = tlo1:thi1;
                    rate1 = 1000*smooth_es(spks1', wid, noedge)';
                    rate_cue(ntr,:) = rate1;
                    [spks2, tlo2, thi2] = spkmtx(spiketrain_temp2(ntr),0,[-1000,3500]);
                    t2 = tlo2:thi2;
                    rate2 = 1000*smooth_es(spks2', wid, noedge)';
                    rate_opp_cue(ntr,:) = rate2;
                end
                for nt = 1:size(rate_cue,2)
                    neuron_ROC(n,nt) = arrayROC(rate_cue(:,nt),rate_opp_cue(:,nt));
                end
            catch
                disp(['error processing neuron  ', num2str(Neurons(n)) '  Dir1=' num2str(Best_class(n))])
            end
        end
        NN1 = size(neuron_ROC,1);
        mr1 = mean(neuron_ROC, 1); % mean rate
        semr1 = std(neuron_ROC, 1)/sqrt(NN1); % SEM
        plt_save.t{g,s} = t1;
        plt_save.Nt{g,s} = NN1;
        plt_save.mr{g,s} = mr1;
        plt_save.semr{g,s} = semr1;
        subplot(2,1,s)
        edges = floor(linspace(1,size(neuron_ROC,1),200))';
        level_t = [];
        for t = 1:size(neuron_ROC,2)
            ROC_t = sort(neuron_ROC(:,t));
            level_t(t,:) = ROC_t(edges);
        end
        level_t = round(level_t,2);
        hm = heatmap(level_t','Colormap',jet,'CellLabelColor','none','ColorLimits',[0.4 1])
        grid off
        Ydata = 1:-0.005:0.001;
        YLabels = string(Ydata);
        YLabels(mod(Ydata,0.2) ~= 0) = " ";
        hm.YDisplayLabels = YLabels;
        Xdata = -1:0.001:3.5;
        XLabels = string(Xdata);
        XLabels(mod(Xdata,0.5) ~= 0) = " ";
        hm.XDisplayLabels = XLabels;
        title(stage{s})
        hm.XDisplayLabels = nan(size(hm.XDisplayData));
        hm.YDisplayLabels = nan(size(hm.YDisplayData));
    end
end