function [ROCarea_all_High, ntrsHigh, err_ntrsHigh,ROCarea_all_Low, ntrsLow, err_ntrsLow] = Neuron_Data_ROCarray_onefileErrtria_MNM(filename,TI,TW,startt,endt,numer,r1r2)

[exam, Classes] = exam_finder_2_classes( filename,r1r2);
classIndex= Classes; %[2 4 10 12]; %adjacent classes
examT = exam; %([1 2 5 6]);% [11. 22.5 45]

fn = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', filename]);
fn_err = [filename '_err'];
fnerr = load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);

time = startt:TI:endt; % xq  %cuetime ==> all time points (1,125)

ROCarea_all_High = nan * ones(1,length(time)); %NaN value matrix
ROCarea_all_Low = nan * ones(1,length(time));
ntrsHigh =nan;
err_ntrsHigh =nan;
ntrsLow = nan;
err_ntrsLow =nan;


%%% Find the max response class %%%
highClass = classIndex(find(examT>0));
lowClass = classIndex(find(examT<0));

%%%%% Find error trials within receptive field (same location as max_class)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntrErrHigh=0;
try
    if ~isempty(highClass)
        for d = 1:length(highClass)
            ntrErrHigh= ntrErrHigh+length(fnerr.MatData.class(highClass(d)).ntr); %extracts the number of trials for that class
        end
    end
catch
end

ntrErrLow=0;
try
    if ~isempty(lowClass)    
        for d = 1:length(lowClass)
            ntrErrLow= ntrErrLow+length(fnerr.MatData.class(lowClass(d)).ntr); %extracts the number of trials for that class
        end
    end
catch
end


if ntrErrHigh > 0 || ntrErrLow>0 % process below only if there are error trials within high Class
    %%%%% Find TS in error trials  %%%%%
    %%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Extract cue rate for the correct trials (max_class_corr) %%%
    n=0;
    for t =(startt-TW/2)/TI:1:(endt-TW/2)/TI % t = (-1-TW/2)/TI:1:(1.5-TW/2)/TI
        n=n+1;
        if  ntrErrHigh >0
            maxcl_cueall = [];
            ntr=[];
            try
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
            catch
            end
            ntrsHigh=sum(ntr);
            
            %%% Extract cue rate for the error trials within receptive field (same location as max_class) %%%
            err_all = [];
            nt=0;
            try
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
            catch
            end
            err_ntrsHigh=nt;
            %         err_all = [err_releaseall_3];%[maxcl_cueall_4 err_releaseall_4];
            if err_ntrsHigh> numer
                ROCarea_all_High(n) = arrayROC(maxcl_cueall, err_all);
            else
                ROCarea_all_High= nan * ones(1,length(time));
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
            try
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
            catch
            end
            ntrsLow=sum(ntr);
            
            %%% Extract cue rate for the error trials within receptive field (same location as max_class) %%%
            err_all = [];
            nt=0;
            try
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
            catch
            end
            err_ntrsLow=nt;
            %         err_all = [err_releaseall_3];%[maxcl_cueall_4 err_releaseall_4];
            if err_ntrsLow> numer
                ROCarea_all_Low(n) = arrayROC(maxcl_cueall, err_all);
            else
                ROCarea_all_Low= nan * ones(1,length(time));
                ntrsLow=NaN;
                err_ntrsLow=NaN;
            end
        end

    end
    
else
    ROCarea_all_High=NaN* ones(1,length(time));
    ntrsHigh=NaN;
    err_ntrsHigh=NaN;
    ROCarea_all_Low=NaN* ones(1,length(time));
    ntrsLow=NaN;
    err_ntrsLow=NaN;
    disp(['No error trials in RF ', filename])
end





