clc; clear all; close all;

NS_BS='NS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  PFC all neurons  %%%%%%%%%%%%%%%%%%%%

[~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_PFC_', NS_BS]);
% [~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['MSNG_PFC_', NS_BS]);
if strcmp(raw1{1,1}, 'Filename')
    raw1=raw1(2:end, :);
end
%%
neuron = 0; 
trial_start = -1000; % in ms
trial_end = 4500; % in ms
spk_for_fano = struct;
for n = 1:length(raw1)
    fn_corr = [raw1{n,1}, '.mat'];
    try
        MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
        disp(fn_corr)
    catch
        try
            fn_corr=[fn_corr(1:9), '1', fn_corr(10:end)];
            MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
            disp(fn_corr)
        catch
            disp('Neuron file not found')
        end
    end
    [max_class_corr, max_class_corr_cuerate, fixrate, max_class_corr_cdrate, distant_class_corr] = MSNG_best_distant(fn_corr);

    if rem(max_class_corr,2)==0
        max_classes=[max_class_corr-1, max_class_corr];
    elseif rem(max_class_corr,2)~=0 && max_class_corr~=17
        max_classes=[max_class_corr, max_class_corr+1];
    elseif max_class_corr==17
        max_classes=max_class_corr;
    end

    if isempty(max_class_corr)
        disp('No Best cue found')
        continue
    end
    
%     Best_class = max_class_corr; 
    Best_class = max_classes;
    if ~isempty(MatData_corr.MatData)
        neuron_data = MatData_corr.MatData.class;
        n_cond = 1;
        for cl = Best_class%1:length(neuron_data)
            try
                spkdata_temp = [];
                [spiketrain_temp, ntrs_temp1] = Get_spiketrain_partial_aligncue(neuron_data(cl).ntr,cl,[trial_start,trial_end]);
                [spkdata_temp, tlo, thi] = spkmtx(spiketrain_temp,0,[trial_start,trial_end]);
                spk_for_fano(n).group(n_cond).spikes = logical(spkdata_temp);
                n_cond = n_cond + 1;
            catch
                disp(['error processing neuron  ', n '  Dir1=' num2str(cl)])
            end
        end
    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end
Result = get_ff(spk_for_fano);
%%
for n=1:length(Result)
    PFC_NS_ff(n,1) = mean(Result(n).FanoFactorAll);
end
%%
clearvars -except PFC_NS_ff NS_BS
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  PPC all neurons  %%%%%%%%%%%%%%%%%%%%

[~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_PPC_', NS_BS]);
% [~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['MSNG_PPC_', NS_BS]);
if strcmp(raw1{1,1}, 'Filename')
    raw1=raw1(2:end, :);
end
%%
neuron = 0; 
trial_start = -1000; % in ms
trial_end = 4500; % in ms
spk_for_fano = struct;
for n = 1:length(raw1)
    fn_corr = [raw1{n,1}, '.mat'];
    try
        MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
        disp(fn_corr)
    catch
        try
            fn_corr=[fn_corr(1:9), '1', fn_corr(10:end)];
            MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
            disp(fn_corr)
        catch
            disp('Neuron file not found')
        end
    end
    [max_class_corr, max_class_corr_cuerate, fixrate, max_class_corr_cdrate, distant_class_corr] = MSNG_best_distant(fn_corr);

    if rem(max_class_corr,2)==0
        max_classes=[max_class_corr-1, max_class_corr];
    elseif rem(max_class_corr,2)~=0 && max_class_corr~=17
        max_classes=[max_class_corr, max_class_corr+1];
    elseif max_class_corr==17
        max_classes=max_class_corr;
    end

    if isempty(max_class_corr)
        disp('No Best cue found')
        continue
    end
    
%     Best_class = max_class_corr; 
    Best_class = max_classes;
    if ~isempty(MatData_corr.MatData)
        neuron_data = MatData_corr.MatData.class;
        n_cond = 1;
        for cl = Best_class%1:length(neuron_data)
            try
                spkdata_temp = [];
                [spiketrain_temp, ntrs_temp1] = Get_spiketrain_partial_aligncue(neuron_data(cl).ntr,cl,[trial_start,trial_end]);
                [spkdata_temp, tlo, thi] = spkmtx(spiketrain_temp,0,[trial_start,trial_end]);
                spk_for_fano(n).group(n_cond).spikes = logical(spkdata_temp);
                n_cond = n_cond + 1;
            catch
                disp(['error processing neuron  ', n '  Dir1=' num2str(cl)])
            end
        end
    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end
Result = get_ff(spk_for_fano);
%%
for n=1:length(Result)
    PPC_NS_ff(n,1) = mean(Result(n).FanoFactorAll);
end

%%
clearvars -except PFC_NS_ff PPC_NS_ff
%%
NS_BS = 'BS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  PFC all neurons  %%%%%%%%%%%%%%%%%%%%

[~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_PFC_', NS_BS]);
% [~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['MSNG_PFC_', NS_BS]);
if strcmp(raw1{1,1}, 'Filename')
    raw1=raw1(2:end, :);
end
%%
neuron = 0; 
trial_start = -1000; % in ms
trial_end = 4500; % in ms
spk_for_fano = struct;
for n = 1:length(raw1)
    fn_corr = [raw1{n,1}, '.mat'];
    try
        MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
        disp(fn_corr)
    catch
        try
            fn_corr=[fn_corr(1:9), '1', fn_corr(10:end)];
            MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
            disp(fn_corr)
        catch
            disp('Neuron file not found')
        end
    end
    [max_class_corr, max_class_corr_cuerate, fixrate, max_class_corr_cdrate, distant_class_corr] = MSNG_best_distant(fn_corr);

    if rem(max_class_corr,2)==0
        max_classes=[max_class_corr-1, max_class_corr];
    elseif rem(max_class_corr,2)~=0 && max_class_corr~=17
        max_classes=[max_class_corr, max_class_corr+1];
    elseif max_class_corr==17
        max_classes=max_class_corr;
    end

    if isempty(max_class_corr)
        disp('No Best cue found')
        continue
    end
    
%     Best_class = max_class_corr; 
    Best_class = max_classes;
    if ~isempty(MatData_corr.MatData)
        neuron_data = MatData_corr.MatData.class;
        n_cond = 1;
        for cl = Best_class%1:length(neuron_data)
            try
                spkdata_temp = [];
                [spiketrain_temp, ntrs_temp1] = Get_spiketrain_partial_aligncue(neuron_data(cl).ntr,cl,[trial_start,trial_end]);
                [spkdata_temp, tlo, thi] = spkmtx(spiketrain_temp,0,[trial_start,trial_end]);
                spk_for_fano(n).group(n_cond).spikes = logical(spkdata_temp);
                n_cond = n_cond + 1;
            catch
                disp(['error processing neuron  ', n '  Dir1=' num2str(cl)])
            end
        end
    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end
Result = get_ff(spk_for_fano);
%%
for n=1:length(Result)
    PFC_BS_ff(n,1) = mean(Result(n).FanoFactorAll);
end
%%
 clearvars -except PFC_BS_ff PFC_NS_ff PPC_NS_ff NS_BS

 %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  PPC all neurons  %%%%%%%%%%%%%%%%%%%%

[~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_PPC_', NS_BS]);
% [~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['MSNG_PPC_', NS_BS]);
if strcmp(raw1{1,1}, 'Filename')
    raw1=raw1(2:end, :);
end
%%
neuron = 0; 
trial_start = -1000; % in ms
trial_end = 4500; % in ms
spk_for_fano = struct;
for n = 1:length(raw1)
    fn_corr = [raw1{n,1}, '.mat'];
    try
        MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
        disp(fn_corr)
    catch
        try
            fn_corr=[fn_corr(1:9), '1', fn_corr(10:end)];
            MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
            disp(fn_corr)
        catch
            disp('Neuron file not found')
        end
    end
    [max_class_corr, max_class_corr_cuerate, fixrate, max_class_corr_cdrate, distant_class_corr] = MSNG_best_distant(fn_corr);

    if rem(max_class_corr,2)==0
        max_classes=[max_class_corr-1, max_class_corr];
    elseif rem(max_class_corr,2)~=0 && max_class_corr~=17
        max_classes=[max_class_corr, max_class_corr+1];
    elseif max_class_corr==17
        max_classes=max_class_corr;
    end

    if isempty(max_class_corr)
        disp('No Best cue found')
        continue
    end
    
%     Best_class = max_class_corr; 
    Best_class = max_classes;
    if ~isempty(MatData_corr.MatData)
        neuron_data = MatData_corr.MatData.class;
        n_cond = 1;
        for cl = Best_class%1:length(neuron_data)
            try
                spkdata_temp = [];
                [spiketrain_temp, ntrs_temp1] = Get_spiketrain_partial_aligncue(neuron_data(cl).ntr,cl,[trial_start,trial_end]);
                [spkdata_temp, tlo, thi] = spkmtx(spiketrain_temp,0,[trial_start,trial_end]);
                spk_for_fano(n).group(n_cond).spikes = logical(spkdata_temp);
                n_cond = n_cond + 1;
            catch
                disp(['error processing neuron  ', n '  Dir1=' num2str(cl)])
            end
        end
    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end
Result = get_ff(spk_for_fano);
%%
for n=1:length(Result)
    PPC_BS_ff(n,1) = mean(Result(n).FanoFactorAll);
end
%%
clearvars -except PPC_BS_ff PPC_NS_ff PFC_BS_ff PFC_NS_ff

%%
p = polyfit(PFC_NS_ff(:,2), PFC_NS_ff(:,1), 1); 
PFC_NS_fano = p(1); 
%%
p = polyfit(PFC_BS_ff(:,2), PFC_BS_ff(:,1), 1); 
PFC_BS_fano = p(1); 
%%
p = polyfit(PPC_NS_ff(:,2), PPC_NS_ff(:,1), 1); 
PPC_NS_fano = p(1); 
%%
p = polyfit(PPC_BS_ff(:,2), PPC_BS_ff(:,1), 1); 
PPC_BS_fano = p(1); 


function Result = get_ff(spk_for_fano)
    %%% COMPUTE FANO FACTOR ETC
    times = 1500:100:4500;
    fanoP.boxWidth = 100; % width of the sliding window in which the counts are made
    fanoP.matchReps = 0;  % number of random choices regarding which points to throw away when matching distributions
    fanoP.binSpacing = 0.25;% 0.25; % bin width when computing distributions of mean counts 0.25
    fanoP.alignTime = 1000; % time of event that data are aligned to (in the output structure, times will be expressed relative to this)
    for g = 1:length(spk_for_fano)
        try
            Result(g) = VarVsMean_XQ(spk_for_fano(g).group, times, fanoP);
        catch
        end
    end
end








