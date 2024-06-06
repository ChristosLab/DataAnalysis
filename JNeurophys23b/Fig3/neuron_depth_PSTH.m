% Depth/layer of PFC neurons in ODR task
% J Zhu, 20230705
%% load data
clearvars
load('sig_odr_data_depth_20230906_raw_max_sulcus.mat');
% load('sig_odr_data_depth_20230608.mat');
file_name = string(neuron_info.Filename);
neuron_info.ID = extract(file_name,1);
%% seg data/label groups
neuron_info.group(neuron_info.Depth<=800)=1;
neuron_info.group(neuron_info.Depth>800&neuron_info.Depth<=1200)=2;
neuron_info.group(neuron_info.Depth>1200)=3;
odr_data = odr_data(~isnan(neuron_info.Depth),:);
neuron_info = neuron_info(~isnan(neuron_info.Depth),:);
% odr_data = odr_data(contains(neuron_info.Stage,'Y'),:);
% neuron_info = neuron_info(contains(neuron_info.Stage,'Y'),:);
%% select neurons
odr_data = odr_data(neuron_info.sulcus==0,:);
neuron_info = neuron_info(neuron_info.sulcus==0,:);
%% PSTH using chronux
trial_start = -1;
trial_end = 3.5;
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
opp_index = [5 6 7 8 1 2 3 4 9];
for g = 1:size(group,1)
    psth_temp = [];
    t_temp = [];
    E_temp = [];
    nn = [];
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_class = neuron_info.cue(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        for n = 1:length(Best_class)
            Opp_class(n) = opp_index(Best_class(n));
        end
        for n = 1:length(Neurons)
            cell_data=odr_data_group(n,:);
            try
                for cl = Best_class(n)
                    temp_data = cell_data{cl};
                    cueon = [temp_data.Cue_onT];    % cueon time for all trials in the class
                    spiketimes = {temp_data.TS};    % spike times for all trials in the class
                    for i = 1:length(cueon)
                        temp_TS = spiketimes{i};
                        temp_TS = temp_TS - cueon(i);
                        temp_spiketimes{i} = temp_TS;    % cueone aligned spike times for all trials in the struct
                    end
                    spikeforchronux = cell2struct(temp_spiketimes,'spiketime',1);
                    [psth_temp(:,n),t_temp(:,n),E_temp(:,n)] = psth(spikeforchronux,0.05,'n',[trial_start,trial_end]);
                end
            catch
                [psth_temp(:,n),t_temp(:,n),E_temp(:,n)] = deal(nan);
            end
        end
        try
            psthmean = mean(psth_temp',1);
            tmean = mean(t_temp',1);
            Emean = std(psth_temp', 1)/sqrt(length(Neurons)); % SEM
        catch
        end
        plt_save.t{g,s} = tmean;
        plt_save.mr{g,s} = psthmean;
        plt_save.semr{g,s} = Emean;
    end
end
disp('finished running')
%% plot
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
        plot([0 0], [0 35], "Color",'k')
        plot([2 2], [0 35], "Color",'k')
    end
    axis([-1 3.5 0 35])
    xlabel('Time (ms)')
    ylabel('Firing rate (spikes/s)')
    title('Spike density \pm 1 SE')
    legend({'adult','','','','young',})
end
origUnits = fig.Units;
fig.Units = fig.PaperUnits;
fig.PaperSize = fig.Position(3:4);% set the Page Size (figure's PaperSize) to match the figure size in the Paper units
fig.Units = origUnits;% restore the original figure Units
%% additional: PSTH on trials using spike density
group = unique(neuron_info.group);
stage = unique(neuron_info.Stage);
plt_save = table;
opp_index = [5 6 7 8 1 2 3 4 9];
for g = 1:size(group,1)
    nn = [];
    for s = 1:size(stage,1)
        Neurons = neuron_info.Neuron(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        Best_class = neuron_info.del(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)));
        odr_data_group = odr_data(neuron_info.group==group(g)&contains(neuron_info.Stage,stage(s)),:);
        for n = 1:length(Best_class)
            Opp_class(n) = opp_index(Best_class(n));
        end
        cdata1 = {};
        ntrs1 = 0;
        nname1 = {};
        for n = 1:length(Neurons)
            try
                [spiketrain_temp1, ntrs_temp1] = Get_spiketrain_partial_aligncue(odr_data_group(n,:),Best_class(n),[]);
                cdata1 = [cdata1; spiketrain_temp1];
                ntrs1 = ntrs1 + ntrs_temp1;
                nname_temp = cell(size(spiketrain_temp1));
                nname_temp(:) = {Neurons(n)};
                nname1 = [nname1; nname_temp];
            catch
                disp(['error processing neuron  ', Neurons(n) '  Dir1=' num2str(Best_class(n))])
            end
        end
        [spks1, tlo, thi] = spkmtx(cdata1,0,[-1000,3500]);
        t1 = tlo:thi;
        t1 = t1/1000;
        wid = 50; % ms
        noedge = 1;
        rate1 = 1000*smooth_es(spks1', wid, noedge)'; % convert to firing rate (spikes/s) by convolving with a Gaussian
        nn1 = size(rate1,1);
        mr1 = mean(rate1, 1); % mean rate
        semr1 = std(rate1, 1)/sqrt(nn1); % SEM
        %             mr_norm = mr/max(mr);
        %             semr_norm = semr/max(mr);
        plt_save.t{g,s} = t1;
        plt_save.Nt{g,s} = nn1;
        plt_save.mr{g,s} = mr1;
        plt_save.semr{g,s} = semr1;
    end
end