clear all; close all; clc;

sheetlist='ODRdistVar'; % ODRdist ODRdistVar MSNG
ExcelFileName=['Database.xlsx'];
pathname=[];
[Neurons_num, Neurons_txt, raw] = xlsread([pathname,ExcelFileName], sheetlist);

warning off MATLAB:divideByZero
if strcmp(Neurons_txt{1,3},'Channel')
    Neurons_txt = Neurons_txt(2:end,:);% original data without title.
    raw = raw(2:end,:);% original data without title.
end
%%
ns=0; bs=0; missing=0; allWidth=[]; neurons=0; no_wf=0;
p1=[]; p2=[]; p3=[];
MSNG=[];
for n =1:length(Neurons_num)
    filename = [Neurons_txt{n,1}];
    fn = [Neurons_txt{n,1},'_',num2str(Neurons_num(n,1))];
    
    load(['R1R2Waveforms\' fn]); %for R1R2 task

%     try  %for MSNG task
%         load(['MSNGWaveforms\' fn]);
%     catch
%         try
%             disp(fn)
%             fn=[fn(1:9) fn(11:end)];
%             load(['MSNGWaveforms\' fn]);
%         catch
%             missing=missing+1;
%             continue
%         end
%     end

    all_waveforms=[];
    for i=1:length(WaveformData.class)
        try
            for j=1:length(WaveformData.class(i).ntr)
                try
                    all_waveforms=[all_waveforms WaveformData.class(i).ntr(j).Waveform];
                catch
                    continue
                end
            end
        catch
            continue
        end
    end
    
    if isempty(all_waveforms)  % checking if there is no waveform information
        no_wf=no_wf+1;
        continue
    end

    mean_wf=mean(all_waveforms,2);
    % normalization
    maxo = max(mean_wf);
    mino = min(mean_wf);
    mean_wf = mean_wf / (maxo-mino);
    mean_wf = mean_wf/max(abs(mean_wf));
    neurons=neurons+1;
    MSNG(neurons).Neurons=fn;
    MSNG(neurons).Channel=raw{n,3};
    MSNG(neurons).Clusters=raw{n,4};
    if strcmpi(sheetlist, 'ODRdist')
        MSNG(neurons).Area=raw{n,6};
    else
        MSNG(neurons).Area=raw{n,5};
    end
%     MSNG(neurons).Waveforms=mean_wf';
    all_avg_wf(neurons,:)=mean_wf';
    mean_wf = mean_wf';
    clear WaveformData all_waveforms
    fprintf('no_wf=%d,   wf=%d, missing=%d, total=%d,  neuron=%d\n',no_wf, neurons, missing, no_wf+neurons+missing, n)
%     save(['Average_waveforms\', fn, '_avg_wf.mat'], 'mean_wf')
end

%%
% load('MSNG_widths_v1.mat')
% load('ODRdist_widths_v1.mat')
load('ODRdistVar_widths_v1.mat')

%%
for i=1: length(MSNG)
    MSNG(i).BinnedWidth=x.Spikewidths_fit(i);
    MSNG(i).SpikeWidth=x.spikewidths_spline(i);

    if MSNG(i).SpikeWidth>300
        MSNG(i).CellType = 'BS';
    else
        MSNG(i).CellType = 'NS';
    end
end
%%
MSNGsigTable = struct2table(MSNG);
xlswrite('NS_BS_Database_20231109.xlsx', table2cell(MSNGsigTable), sheetlist);





