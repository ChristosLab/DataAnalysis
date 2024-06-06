clc; clear all; close all;

[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109.xlsx',['AllNeurons']);

%%
neuron = 0;
for n=1:length(raw)
    n
    filename             = raw{n, 1}
    wf_duration          = raw{n, 6};
%     if wf_duration == 725
%         continue
%     end
    [fixrates, fix_rate] = baseline_firingrate(filename);
    if ~isempty(fix_rate)
        neuron = neuron + 1;
        all_wf(neuron, :) = [wf_duration, fix_rate];
        all_fixrates(neuron, 1) = {fixrates};
    end

end

%% all datasets
fig1 = figure;
widths = all_wf(:,1);
fix_rates = all_wf(:,2);
ns_i = find(widths<=300);
bs_i = find(widths>300);
[h, p, ci, stats] = ttest2(fix_rates(ns_i), fix_rates(bs_i))

% [h, p, ci, stats] = ttest2(group1, group2);
for n = 1:length(all_wf)
    if all_wf(n,1) >300
        plot(all_wf(n,1), log10(all_wf(n,2)), 'Marker', 'o', 'Color', 'b', 'MarkerSize', 5); hold on; 
%         plot(all_wf(n,1), all_wf(n,2), 'Marker', '.', 'Color', 'b', 'MarkerSize', 10); hold on; 
    else
        plot(all_wf(n,1), log10(all_wf(n,2)), 'Marker', 'o', 'Color', 'r', 'MarkerSize', 5); hold on;
%         plot(all_wf(n,1), all_wf(n,2), 'Marker', '.', 'Color', 'r', 'MarkerSize', 10); hold on;
    end
end
% % linear regression
% y = log10(all_wf(:,2));
% x = all_wf(:,1);
% y = y(find(~isinf(y)));
% x = x(find(~isinf(y)));
% 
% ns_i = find(x<=300);
% bs_i = find(x>300);
% [h, p, ci, stats] = ttest2(y(ns_i), y(bs_i))
% 
% coefficients = polyfit(x,  y, 1);
% x_fit = linspace(min(x), max(x), 1000);
% y_fit = polyval(coefficients, x_fit);
% 
% % Plot the linear regression line
% plot(x_fit, y_fit, 'k', 'LineWidth', 2);
% add labels
xlabel('Waveform Duration ($\mu s$)', 'FontSize', 15, 'Interpreter','latex');
ylabel('Log fixation firing rate (sp/s)', 'FontSize', 15, 'Interpreter','latex');
title('MSNG + ODRdist + ODRdistVar', 'FontSize', 15, 'Interpreter','latex')

set(fig1, 'Renderer', 'painters');
exportgraphics(fig1, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure3\scatterplot_fixrateVSwidth.emf', 'ContentType', 'vector');
print(fig1, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure3\scatterplot_fixrateVSwidth.pdf', '-vector', '-bestfit', '-dwinc');


%% histogram plot
duros=all_wf(:,1)';
NS =sum(duros<=300);
BS=sum(duros>300);
fig2=figure; 
hist=histogram(duros',[min(min(duros),150):25:max(duros)], 'Normalization', 'count','EdgeColor', 'k', 'FaceColor', 'blue');
% title(sprintf('Classification with normalized detrended amplitude (NS=%d, BS=%d, Unknown=%d)', ns,bs,unknown))
x = 300;  % Replace with your threshold value

% Set different colors for bars above and below the threshold
hold on;
binWidth = hist.BinWidth;
bar_centers = hist.BinEdges(1:end-1) + 12.5;
bar_heights = hist.Values;
above_threshold = bar_centers >= x;

% colors={[0.6 0.6 1] , [1 0.6 0.6]} ;
colors={'b' , 'r'} ;
b1=bar(bar_centers(above_threshold), bar_heights(above_threshold), 'BarWidth', 1,'FaceColor', colors{1}, 'EdgeColor', 'k');
b2=bar(bar_centers(~above_threshold), bar_heights(~above_threshold), 'BarWidth', 1, 'FaceColor', colors{2}, 'EdgeColor', 'k');



hold off;

legend([b1, b2], 'Broad units', 'Narrow units', 'Interpreter', 'latex', 'FontSize',15);  % Show legend
set(gca, 'FontSize', 15);
xlabel('Waveform Duration ($\mu$s)', 'Interpreter', 'latex', 'FontSize',15); 
ylabel('Number of Neurons', 'Interpreter', 'latex', 'FontSize',15);
title(sprintf('Histogram of All Widths (NS=%d, BS=%d)', NS, BS), 'Interpreter', 'latex', 'FontSize',15)
box off
%%
duros=duros(~isnan(duros));
[dip,p,xlow,xup] = HartigansDipSignifTest(duros,10000);
fprintf('dip value: %f    p-value: %f\n ',dip,p);

%%
set(fig2, 'Renderer', 'painters');s
exportgraphics(fig2, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure3\NS_BS_histogram_20231109.emf', 'ContentType', 'vector');
print(fig2, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure3\NS_BS_histogram_20231109.pdf', '-vector', '-bestfit', '-dwinc');
%% MSNG
figure
for n = 1:1436
    if all_wf(n,1) >300
        plot(all_wf(n,1), all_wf(n,2), 'Marker', 'square', 'Color', 'b', 'MarkerSize', 10); hold on; 
    else
        plot(all_wf(n,1), all_wf(n,2), 'Marker', 'v', 'Color', 'r', 'MarkerSize', 10); hold on;
    end
end
xlabel('Waveform Duration ($\mu s$)', 'FontSize', 15, 'Interpreter','latex');
ylabel('Log fixation firing rate (sp/s)', 'FontSize', 15, 'Interpreter','latex');
title('MSNG', 'FontSize', 15, 'Interpreter','latex')
%% ODRdistVar
figure
for n = 1437:1880
    if all_wf(n,1) >300
        plot(all_wf(n,1), all_wf(n,2), 'Marker', 'square', 'Color', 'b', 'MarkerSize', 10); hold on; 
    else
        plot(all_wf(n,1), all_wf(n,2), 'Marker', 'v', 'Color', 'r', 'MarkerSize', 10); hold on;
    end
end
xlabel('Waveform Duration ($\mu s$)', 'FontSize', 15, 'Interpreter','latex');
ylabel('Fixation firing rate (sp/s)', 'FontSize', 15, 'Interpreter','latex');
title('ODRdistVar', 'FontSize', 15, 'Interpreter','latex')
%% ODRdist
figure
for n = 1881:2461
    if all_wf(n,1) >300
        plot(all_wf(n,1), all_wf(n,2), 'Marker', 'square', 'Color', 'b', 'MarkerSize', 10); hold on; 
    else
        plot(all_wf(n,1), all_wf(n,2), 'Marker', 'v', 'Color', 'r', 'MarkerSize', 10); hold on;
    end
end
xlabel('Waveform Duration ($\mu s$)', 'FontSize', 15, 'Interpreter','latex');
ylabel('Fixation firing rate (sp/s)', 'FontSize', 15, 'Interpreter','latex');
title('ODRdist', 'FontSize', 15, 'Interpreter','latex')
%%
function [all_fixrates, fixrate]    = baseline_firingrate(filename)

    fixrate      = [];
    all_fixrates = [];
    dir          = 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\';

    try
        fn = [dir, filename];
        load(fn)
    catch
        try
            fn = [dir, filename(1:9), '1', filename(10:end)];
            load(fn)
        catch
            disp('Neuron file not found')
            return
        end
    end

    if ~isempty(MatData)
        for c = 1:length(MatData.class)
            if ~isempty(MatData.class(c).ntr)
                all_fixrates = [all_fixrates, mean([MatData.class(c).ntr.fix])];
            else
                disp('Empty class!!')
            end
        end
        fixrate = mean(all_fixrates);
    else
        disp('Empty Matdata')
    end
end



















