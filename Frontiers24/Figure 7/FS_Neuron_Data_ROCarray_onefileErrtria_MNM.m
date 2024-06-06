function [ROCarea_all_High, ntrsHigh, err_ntrsHigh,ROCarea_all_Low, ntrsLow, err_ntrsLow] = Neuron_Data_ROCarray_onefileErrtria_MNM(fname,TI,TW,startt,endt,numer)
% filename,TI,TW,cWindow(1),cWindow(2),numer,phase,bigRF
% function [ROCarea_all_1 ROCarea_all_2 ROCarea_all_3 ROCarea_all_4] = Neuron_Data_ROCarray_onefileDifftaskErrortrial(filename,TI,TW)
% Time-resolved ROC analysis from Saliency experiment.
% Trials from arrays of different colors are pooled together.
% TI, TW parameters define size and overlap of sliding windows.
%
% This is for the diffcluty tasks (e.g. NoJump32_4DC_EW_MNM1rec) with 2
% locations (EW/NS), 2 cue colors (green/red), 4 difficulty levels and
% 2 task types (Match/NonMatch1).
% Edited: Apr 29 2010
%
% Modified for error analysis.
% Use only Level 3 trials in diffictul task. Load neurondata file for both correct and
% error trials.
% Find RF (max response location) using level 1 trials and collect TS for correct (got location) and error (missed location)
% trials for Level 3 in the receptive field reagardless of type of error (during match/nonmatch classes, hold/release errors).
% Run ROC analysis on correct vs. error firing rate.
% Modified to use neurons with more than n (>n) error trials.
% Modified to switch case for 'all' (using both match and nm1 trials),'match'(collecting only immediate match
% trials) or 'nm1' (collecting only nm1 trials).
% Last edited: Jul 3 2013

% TI = 0.01; % time interval(sec) for ROC
% TW = 0.05; % time window to collect firing rates
% cuetime =  -0.25:TI:0.75; %-1:TI:1.5;
cuetime = startt:TI:endt; % xq
crit=0.1;
% [ exam ] = Neuron_Data_average_CorrectVSIncorrectFR_loc17_crit_exam( filename, crit );%[2:2:16]
[ exam ] = FS_Neuron_Data_average_CorrectVSIncorrectFR_loc17_crit( fname, crit );
% classIndexT=[2:2:16];
classIndex= [2 4 10 12];%[2 10];
if ~isempty(exam)
    examT = exam([1 2 5 6]);% exam([1 5]);% [11.25 22.5 -11.25 -22.5]
else
    ROCarea_all_High=NaN* ones(1,length(cuetime));
    ntrsHigh=NaN;
    err_ntrsHigh=NaN;
    ROCarea_all_Low=NaN* ones(1,length(cuetime));
    ntrsLow=NaN;xticklabels({'A', 'B', 'C', 'D', 'E', 'F', 'G'})
    err_ntrsLow=NaN;
    return
end



% cd 'C:\Users\chungs6\Documents\MATLAB\Independent Data Analysis Project\replicate li et al analyses\allMNMcorrect'
try
    fn = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', fname]);
catch
    fn = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', [fname(1:9),'1', fname(10:end)]]);
end
% cd 'C:\Users\chungs6\Documents\MATLAB\Independent Data Analysis Project\replicate li et al analyses\allMNMerror'
try
    fnerr = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', [fname 'err']]);
catch
    try
        fnerr = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', [fname(1:9),'1', fname(10:end), 'err']]);
    catch
        disp([fname 'err not found'])
    end
end
% fnbhv = load([filename(1:8)]);

ROCarea_all_High = nan * ones(1,length(cuetime));
ROCarea_all_Low = nan * ones(1,length(cuetime));
ntrsHigh =nan;
err_ntrsHigh =nan; 
ntrsLow = nan;
err_ntrsLow =nan;



%%% Find the max response class %%%
highClass = classIndex(find(examT>0));
lowClass = classIndex(find(examT<0));

% if ~isempty(fn.MatData)&&~isempty(fnerr.MatData)
%     for b = 1:8%length(fn.MatData.class) % Number of classes in one difficulty level
%         eval(['temprate(b) = mean([fn.MatData.class(b).ntr.' phase ']);']);
%     end
%     max_class = find(temprate == max(temprate)); % Find the class with maximum mean cuerates.
%     if bigRF
%         tt=circshift([1:8],[1 5-max_class]);%length(fn.MatData.class)
%         max_classall = tt(4:6);% pick center 3 classes
%     else
%         max_classall = max_class;
%     end
% else
%     disp(['empty ', filename])
%     ROCarea_all=NaN * ones(1,length(cuetime));
%     ntrs=NaN;
%     err_ntrs=NaN;
%     return
% end
% 
%%%%% Find error trials within receptive field (same location as max_class)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntrErrHigh=0;
if ~isempty(highClass)
    
    for d = 1:length(highClass)
        %         emperr(d) = isempty(fnerr.MatData.class(d).ntr);
        ntrErrHigh= ntrErrHigh+length(fnerr.MatData.class(highClass(d)).ntr);
    end
    %     fnd_errc = find(emperr==0); % all error classes
    %     errcl =intersect(fnd_errc,max_classall); % only error classes within RF
end
ntrErrLow=0;
if ~isempty(lowClass)    
    for d = 1:length(lowClass)
        %         emperr(d) = isempty(fnerr.MatData.class(d).ntr);
        ntrErrLow= ntrErrLow+length(fnerr.MatData.class(lowClass(d)).ntr);
    end
%     fnd_errc = find(emperr==0); % all error classes
%     errcl =intersect(fnd_errc,max_classall); % only error classes within RF
end


if ntrErrHigh > 0 || ntrErrLow>0% process below only if there are error trials within high Class
    %%%%% Find TS in error trials  %%%%%
    %%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Extract cue rate for the correct trials (max_class_corr) %%%
    n=0;
    for t =(startt-TW/2)/TI:1:(endt-TW/2)/TI % t = (-1-TW/2)/TI:1:(1.5-TW/2)/TI
        n=n+1;
        if  ntrErrHigh >0
            maxcl_cueall = [];
            ntr=[];
            for nc = 1:length(highClass)
                maxClass = highClass(nc);
                for nn = 1:length(fn.MatData.class(maxClass).ntr)
                    ntr_TS1 = fn.MatData.class(maxClass).ntr(nn).TS;
                    Cue_onT1 = fn.MatData.class(maxClass).ntr(nn).Cue_onT;
                    maxcl_TS1 = ntr_TS1(find(ntr_TS1 >= Cue_onT1+(TI*t) & ntr_TS1 < Cue_onT1+TI*t+TW));
                    maxcl_cue1 = length(maxcl_TS1);
                    maxcl_cueall = [maxcl_cueall maxcl_cue1];
                end
                ntr=[ntr nn];
            end
            ntrsHigh=sum(ntr);
            
            %%% Extract cue rate for the error trials within receptive field (same location as max_class) %%%
            err_all = [];
            nt=0;
            for ne = 1:length(highClass) % Find TS in error trials within RF
                for ii = 1:length(fnerr.MatData.class(highClass(ne)).ntr)
                    ntr_TS_err1 = fnerr.MatData.class(highClass(ne)).ntr(ii).TS;
                    Cue_onT_err1 = fnerr.MatData.class(highClass(ne)).ntr(ii).Cue_onT;
                    err_TS1 = ntr_TS_err1(find(ntr_TS_err1 >= Cue_onT_err1+(TI*t) & ntr_TS_err1 < Cue_onT_err1+TI*t+TW));
                    err_cue1 = length(err_TS1);
                    err_all = [err_all err_cue1];
                    nt=nt+1;
                end
            end
            err_ntrsHigh=nt;
            %         err_all = [err_releaseall_3];%[maxcl_cueall_4 err_releaseall_4];
%             cd 'C:\Users\chungs6\Documents\MATLAB\Independent Data Analysis Project\ROC_analysis_FS'
            if err_ntrsHigh> numer
                ROCarea_all_High(n) = arrayROC(maxcl_cueall, err_all);
            else
                ROCarea_all_High= nan * ones(1,length(cuetime));
                ntrsHigh=NaN;
                err_ntrsHigh=NaN;
%                 disp(['Insufficient number of error trials  ', filename])
%                 return
            end
        end
        % low calss
        if  ntrErrLow >0
            maxcl_cueall = [];
            ntr=[];
            for nc = 1:length(lowClass)
                maxClass = lowClass(nc);
                for nn = 1:length(fn.MatData.class(maxClass).ntr)
                    ntr_TS1 = fn.MatData.class(maxClass).ntr(nn).TS;
                    Cue_onT1 = fn.MatData.class(maxClass).ntr(nn).Cue_onT;
                    maxcl_TS1 = ntr_TS1(find(ntr_TS1 >= Cue_onT1+(TI*t) & ntr_TS1 < Cue_onT1+TI*t+TW));
                    maxcl_cue1 = length(maxcl_TS1);
                    maxcl_cueall = [maxcl_cueall maxcl_cue1];
                end
                ntr=[ntr nn];
            end
            ntrsLow=sum(ntr);
            
            %%% Extract cue rate for the error trials within receptive field (same location as max_class) %%%
            err_all = [];
            nt=0;
            for ne = 1:length(lowClass) % Find TS in error trials within RF
                for ii = 1:length(fnerr.MatData.class(lowClass(ne)).ntr)
                    ntr_TS_err1 = fnerr.MatData.class(lowClass(ne)).ntr(ii).TS;
                    Cue_onT_err1 = fnerr.MatData.class(lowClass(ne)).ntr(ii).Cue_onT;
                    err_TS1 = ntr_TS_err1(find(ntr_TS_err1 >= Cue_onT_err1+(TI*t) & ntr_TS_err1 < Cue_onT_err1+TI*t+TW));
                    err_cue1 = length(err_TS1);
                    err_all = [err_all err_cue1];
                    nt=nt+1;
                end
            end
            err_ntrsLow=nt;
            %         err_all = [err_releaseall_3];%[maxcl_cueall_4 err_releaseall_4];
%             cd 'C:\Users\chungs6\Documents\MATLAB\Independent Data Analysis Project\ROC_analysis_FS'
            if err_ntrsLow> numer
                ROCarea_all_Low(n) = arrayROC(maxcl_cueall, err_all);
            else
                ROCarea_all_Low= nan * ones(1,length(cuetime));
                ntrsLow=NaN;
                err_ntrsLow=NaN;
%                 disp(['Insufficient number of error trials in RF ', filename])
%                 return
            end
        end

    end
    
else
    ROCarea_all_High=NaN* ones(1,length(cuetime));
    ntrsHigh=NaN;
    err_ntrsHigh=NaN;
     ROCarea_all_Low=NaN* ones(1,length(cuetime));
    ntrsLow=NaN;xticklabels({'A', 'B', 'C', 'D', 'E', 'F', 'G'})
    err_ntrsLow=NaN;
    disp(['No error trials in RF ', fname])
end
% figure
% plot(cuetime,ROCarea_all_3);
% title(['ROC','  ',filename(1:6),'\_',filename(8),'\_',filename(10:length(filename)),'Level 3'])
% print('-djpeg', [filename,'.jpg'])




% close all
