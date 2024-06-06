function [ROCarea, ntrsBest, ntrsDiag] = ROCclass_best_diag(filename,TI,TW,startt,endt,r1r2,CorrErr) 
% function [ROCarea, ntrsBest, ntrsDiag, ROCarea_best, ROCarea_diag] = ROC_best_diag_err(filename,TI,TW,startt,endt,r1r2)

[exam, Classes] = exam_finder_2_classes( filename,r1r2);
classIndex= Classes; 
examT = exam; 

if CorrErr==1
    fn = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', filename]);
elseif CorrErr==0
    fn_err = [filename '_err'];
    fn = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);
end

time = startt:TI:endt; % xq  %cuetime ==> all time points (1,125)

ROCarea = nan*ones(1,length(time)); %NaN value matrix
% ROCarea_best = nan*ones(1,length(time));
% ROCarea_diag = nan*ones(1,length(time));

ntrsBest =nan;
ntrsDiag =nan;

%%% Find the max response class %%%
BestClass = classIndex(find(examT>0));
DiagClass = classIndex(find(examT<0));

%% Find  trials within receptive field and in diagonal location %%
ntrsBest=0;
try
    if ~isempty(BestClass)
        for d = 1:length(BestClass)
            ntrsBest= ntrsBest+length(fn.MatData.class(BestClass(d)).ntr); %extracts the number of trials for that class
        end
    end
catch
end

ntrsDiag=0;
try
    if ~isempty(DiagClass)
        for d = 1:length(DiagClass)
            ntrsDiag= ntrsDiag+length(fn.MatData.class(DiagClass(d)).ntr); %extracts the number of trials for that class
        end
    end
catch
end


if ntrsBest > 0 && ntrsDiag>0 % process below only if there are trials within Best and Diag Classes
    
    %%% Extract cue rate for the correct trials (max_class_corr) %%%
    n=0;
    for t =(startt-TW/2)/TI:1:(endt-TW/2)/TI % t = (-1-TW/2)/TI:1:(1.5-TW/2)/TI
        n=n+1;
        bestcl_cueall = [];
        ntr=[];
        try
            for nc = 1:length(BestClass)
                maxClass = BestClass(nc);
                maxcl_cueall=[];
                for nn = 1:length(fn.MatData.class(maxClass).ntr)
                    ntr_TS1 = fn.MatData.class(maxClass).ntr(nn).TS;
                    Cue_onT1 = fn.MatData.class(maxClass).ntr(nn).Cue_onT;
                    maxcl_TS1 = ntr_TS1(find(ntr_TS1 >= Cue_onT1+(TI*t) & ntr_TS1 < Cue_onT1+TI*t+TW));
                    maxcl_cue1 = length(maxcl_TS1);
                    maxcl_cueall = [maxcl_cueall maxcl_cue1];
%                     bestcl_cueall = [bestcl_cueall maxcl_cue1];
                end
                if ~isempty(maxcl_cueall)
                    bestcl_cueall=[bestcl_cueall mean(maxcl_cueall)];
                end
                ntr=[ntr nn];
            end
        catch
        end
        ntrsHigh=sum(ntr);

        %%
        diagcl_cueall = [];
        ntr=[];
        try
            for nc = 1:length(DiagClass)
                minClass = DiagClass(nc);
                mincl_cueall=[];
                for nn = 1:length(fn.MatData.class(minClass).ntr)
                    ntr_TS1 = fn.MatData.class(minClass).ntr(nn).TS;
                    Cue_onT1 = fn.MatData.class(minClass).ntr(nn).Cue_onT;
                    mincl_TS1 = ntr_TS1(find(ntr_TS1 >= Cue_onT1+(TI*t) & ntr_TS1 < Cue_onT1+TI*t+TW));
                    mincl_cue1 = length(mincl_TS1);
                    mincl_cueall = [mincl_cueall mincl_cue1];
%                     diagcl_cueall = [diagcl_cueall mincl_cue1];
                end
                if ~isempty(mincl_cueall)
                    diagcl_cueall=[diagcl_cueall mean(mincl_cueall)];
                end
                ntr=[ntr nn];
            end
        catch
        end
        ntrsLow=sum(ntr);
        ROCarea(n) = arrayROC(bestcl_cueall, diagcl_cueall);
%         ROCarea_best(n) = arrayROC(bestcl_cueall, 0.5*ones(1, length(bestcl_cueall)));
%         ROCarea_diag(n) = arrayROC(diagcl_cueall, 0.5*ones(1, length(diagcl_cueall)));

    end
else
    ROCarea=NaN* ones(1,length(time));
    ntrsHigh=NaN;
    ntrsDiag=NaN;
    disp(['Not enough trials in RF ', filename])
end





