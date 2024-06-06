% Depth/layer of PFC neurons in ODR task
% J Zhu, 20230606
%% load data
clearvars
load('sig_odr_data_depth_20230829_raw_max_sulcus.mat');
file_name = string(neuron_info.Filename);
neuron_info.ID = extract(file_name,1);
%% select neuron: optional
odr_data = odr_data(neuron_info.sulcus==0,:);
neuron_info = neuron_info(neuron_info.sulcus==0,:);
%% plot destribution of neurons depth
DI_plt1 = -neuron_info.Depth;
group_plt1 = neuron_info.ID;
id = ['i','j','k','l'];
hf = figure;
allTS = [];
for m = 1:length(id)
    plt = DI_plt1(group_plt1==id(m));
    subplot(1,2,1)
    hold on
    TS = plt;
    for o = 1:length(TS)
        line([m-.9,m-.1],[TS(o),TS(o)],'Color', 'k')
    end
    line([m-1 m],[0 0],'Color','b')
    line([m-1 m],[-800 -800],'Color','b')
    line([m-1 m],[-1200 -1200],'Color','b')
    allTS = [allTS; TS];
    axis([0 m -3500 100])
    hold off
end
subplot(1,2,2)
hold on
bin_width = 200;
bin_edges = min(allTS):bin_width:0;
psth = histc(allTS,bin_edges)/size(allTS,1);
bins = bin_edges+0.5*bin_width;
hb = barh(bins,psth);
set(hb,'EdgeColor','k');
set(hb,'FaceColor','k');
axis([0 0.2 -3500 100])
hold off
%% seg data/label groups
neuron_info.group(neuron_info.Depth<=800)=1;
neuron_info.group(neuron_info.Depth>800&neuron_info.Depth<=1200)=2;
neuron_info.group(neuron_info.Depth>1200)=3;
odr_data = odr_data(~isnan(neuron_info.Depth),:);
neuron_info = neuron_info(~isnan(neuron_info.Depth),:);
% odr_data = odr_data(contains(neuron_info.Stage,'Y'),:);
% neuron_info = neuron_info(contains(neuron_info.Stage,'Y'),:);
%% PSTH
group = unique(neuron_info.group);
plt_save = table;
for g = 1:size(group,1)
    Neurons = neuron_info.Neuron(neuron_info.group==group(g));
    Best_class = neuron_info.cue(neuron_info.group==group(g));
    odr_data_group = odr_data(neuron_info.group==group(g),:);
    opp_index = [5 6 7 8 1 2 3 4 9];
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
    plt_save.t{g} = t1;
    plt_save.Nt{g} = nn1;
    plt_save.mr{g} = mr1;
    plt_save.semr{g} = semr1;
end
%% plot
% define a few of colors
my_color = linspecer(3);
my_color = my_color([1:3],:);

figure
clf
set(gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.2,0.2,0.3,0.4] );
for g = 1:size(group,1)
    xlo = min(plt_save.t{g});
    xhi = max(plt_save.t{g});
    ylo = 0;
    yhi = max(plt_save.mr{g})*1.1;

    hold on;
    set(gca, 'tickdir', 'out')

    % plot errors as a shaded area
    y11 = plt_save.mr{g} + plt_save.semr{g};
    y12 = plt_save.mr{g} - plt_save.semr{g};
    fill([plt_save.t{g} plt_save.t{g}(end:-1:1)], [y11 y12(end:-1:1)], 'w', ...
        'EdgeColor', my_color(g,:), 'FaceColor', my_color(g,:),'FaceAlpha',0.3,'EdgeAlpha',0.3)

    % plot mean on top
    plot(plt_save.t{g}, plt_save.mr{g}, 'Color',my_color(g,:), 'linewidth', 1)
    plot([0 0], ylim, "Color",my_color(g,:))
end
axis([-1 3.5 0 35])
xlabel('Time (ms)')
ylabel('Firing rate (spikes/s)')
title('Spike density \pm 1 SE')
legend({'superficial','','','mid','','','deep','',''})
%% additional: PSTH avg on neuron
group = unique(neuron_info.group);
plt_save = table;
for g = 1:size(group,1)
    Neurons = neuron_info.Neuron(neuron_info.group==group(g));
    Best_class = neuron_info.cue(neuron_info.group==group(g));
    odr_data_group = odr_data(neuron_info.group==group(g),:);
    opp_index = [5 6 7 8 1 2 3 4 9];
    for n = 1:length(Best_class)
        Opp_class(n) = opp_index(Best_class(n));
    end
    mr_neuron=[];
    for n = 1:length(Neurons)
        try
            [spiketrain_temp1, ntrs_temp1] = Get_spiketrain_partial_aligncue(odr_data_group(n,:),Best_class(n),[]);
            cdata1 = [cdata1; spiketrain_temp1];
            ntrs1 = ntrs1 + ntrs_temp1;
            nname_temp = cell(size(spiketrain_temp1));
            nname_temp(:) = {Neurons(n)};
            nname1 = [nname1; nname_temp];

            [spks1, tlo, thi] = spkmtx(cdata1,0,[-1000,3500]);
            t1 = tlo:thi;
            wid = 50; % ms
            noedge = 1;
            rate1 = 1000*smooth_es(spks1', wid, noedge)'; % convert to firing rate (spikes/s) by convolving with a Gaussian
            mr_neuron(n,:) = mean(rate1,1);
        catch
            disp(['error processing neuron  ', Neurons(n) '  Dir1=' num2str(Best_class(n))])
        end
    end
    nn1 = size(mr_neuron,1);
    mr1 = mean(mr_neuron, 1); % mean rate
    semr1 = std(mr_neuron, 1)/sqrt(nn1); % SEM
    plt_save.t{g} = t1;
    plt_save.Nt{g} = nn1;
    plt_save.mr{g} = mr1;
    plt_save.semr{g} = semr1;
end
%% additional: PSTH using chronux (sanity check)
trial_start = -1;
trial_end = 3.5;
group = unique(neuron_info.group);
plt_save = table;
for g = 1:size(group,1)
    psth_temp = [];
    t_temp = [];
    E_temp = [];
    Neurons = neuron_info.Neuron(neuron_info.group==group(g));
    Best_class = neuron_info.del(neuron_info.group==group(g));
    odr_data_group = odr_data(neuron_info.group==group(g),:);
    opp_index = [5 6 7 8 1 2 3 4 9];
    for n = 1:length(Best_class)
        Opp_class(n) = opp_index(Best_class(n));
    end
    cdata1 = {};
    ntrs1 = 0;
    nname1 = {};
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
        psthmean = mean(psth_temp,2,"omitnan");
        tmean = mean(t_temp,2,"omitnan");
        Emean = mean(E_temp,2,"omitnan");
    catch
    end
    plt_save.t{g} = tmean';
    plt_save.mr{g} = psthmean';
    plt_save.semr{g} = Emean';
end
disp('finished running')