function [ROCarea, ntrsHigh, ntrsDiag] = ROCntr_best_worst(filename,TI,TW,startt,endt)

try
    fn = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', filename]);
catch
    fn = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', filename(1:9), '1', filename(10:end)]);
    filename = [filename(1:9), '1', filename(10:end)];
end

[exam, Classes] = exam_finder_MSNG( filename);
classIndex= Classes; 
examT = exam; 


time = startt:TI:endt; % xq  %cuetime ==> all time points (1,125)

ROCarea = nan*ones(1,length(time)); %NaN value matrix

ntrsBest =nan;
ntrsDiag =nan;

%%% Find the max response class %%%
BestClass = classIndex(find(examT>0));
DiagClass = classIndex(find(examT<0));

%%%%% Find  trials within receptive field and in diagonal location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                for nn = 1:length(fn.MatData.class(maxClass).ntr)
                    ntr_TS1 = fn.MatData.class(maxClass).ntr(nn).TS;
                    Cue_onT1 = fn.MatData.class(maxClass).ntr(nn).Cue_onT;
                    maxcl_TS1 = ntr_TS1(find(ntr_TS1 >= Cue_onT1+(TI*t) & ntr_TS1 < Cue_onT1+TI*t+TW));
                    maxcl_cue1 = length(maxcl_TS1);
                    bestcl_cueall = [bestcl_cueall maxcl_cue1];
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
                maxClass = DiagClass(nc);
                for nn = 1:length(fn.MatData.class(maxClass).ntr)
                    ntr_TS1 = fn.MatData.class(maxClass).ntr(nn).TS;
                    Cue_onT1 = fn.MatData.class(maxClass).ntr(nn).Cue_onT;
                    maxcl_TS1 = ntr_TS1(find(ntr_TS1 >= Cue_onT1+(TI*t) & ntr_TS1 < Cue_onT1+TI*t+TW));
                    maxcl_cue1 = length(maxcl_TS1);
                    diagcl_cueall = [diagcl_cueall maxcl_cue1];
                end
                ntr=[ntr nn];
            end
        catch
        end
        ntrsLow=sum(ntr);
        
        ROCarea(n) = arrayROC(bestcl_cueall, diagcl_cueall);
%         ROCarea_best(n) = arrayROC(bestcl_cueall, 0.5*ones(1,length(bestcl_cueall)));
%         ROCarea_diag(n) = arrayROC(diagcl_cueall,  0.5*ones(1,length(diagcl_cueall)));

    end
else
    ROCarea=NaN* ones(1,length(time));
    ntrsHigh=NaN;
    ntrsDiag=NaN;
    disp(['Not enough trials in RF ', filename])
end
