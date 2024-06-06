%%Results are in result, frpool, fr_alltrials_err;
%% check excel_output for data result
%% extract delay firing rate from class 17 as 0 degree firing rate;
function [ op ] = Neuron_Data_average_CorrectVSIncorrectFR_loc17_crit( fname, crit )
% filename = 'lem589_2_8410';

try
    load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', fname]);
catch
    load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', [fname(1:9),'1', fname(10:end)]]);
end
time_period = 1;%% test last 0.5/1/3 delay period firing rate
% fr = zeros(1,16);
% location = [-90 -45 -22.5 -11.25 11.25 22.5 45 90]
% crit = 0; %%criterion for firing rate

%%generate delay rate for different locations from NM trials
for z = 2:2:16 %%only +11.25 and -11.25
    try
        if length(MatData.class(z).ntr) == 0
        anmrate = 0;
        else
            for i = 1:length(MatData.class(z).ntr)
                ts = MatData.class(z).ntr(i).TS;
                a = find(ts >  MatData.class(z).ntr(i).Sample_onT-time_period & ts < MatData.class(z).ntr(i).Sample_onT);
                anmrate(i) = length(a)/time_period;
            end
        end
    catch
        anmrate = NaN;
    end
    %     rate(c/2) = ((sum(mrate)+sum(nmrate))/(length(MatData.class(c).ntr)+length(MatData.class(c-1).ntr)))
    sumrate(z/2) = sum(anmrate);
    try
        sumtrialnum(z/2) = length(MatData.class(z).ntr);
        stdrate(z/2) = std(anmrate);
    catch
        sumtrialnum(z/2) = NaN;
        stdrate(z/2) = NaN;
    end
    clear anmrate
end
delay_rate = sumrate./sumtrialnum;

% just for FS analysis 
try 
    for j = 1:length(MatData.class(16).ntr)
        ts17 = MatData.class(16).ntr(j).TS;
        d = find(ts17 >  MatData.class(16).ntr(j).Sample_onT-time_period & ts17 < MatData.class(16).ntr(j).Sample_onT);
        zerorate(j) = length(d)/time_period;
    end
averagefr = mean(zerorate);
catch
    averagefr = NaN;
end

% for j = 1:length(MatData.class(17).ntr)
%     ts17 = MatData.class(17).ntr(j).TS;
%     d = find(ts17 >  MatData.class(17).ntr(j).Sample_onT-time_period & ts17 < MatData.class(17).ntr(j).Sample_onT);
%     zerorate(j) = length(d)/time_period;
% end
% averagefr = mean(zerorate);%% averagefr indicates firing rate from class 17

clear ts a ts17 d

%%generate delay rate for different locations for error trials
% cd 'C:\Users\chungs6\Documents\MATLAB\Independent Data Analysis Project\replicate li et al analyses\allMNMerror'
try
    load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', [fname 'err']]);
catch
    try
        load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', [fname(1:9),'1', fname(10:end), 'err']]);
    catch
        disp([filename 'err not found']);
        exam=[];op = [exam];
        return
    end
end
err_delay_rate = nan(1,8);
erranmrate = [];
try 
    for y = 2:2:length(MatData.class) %%only +11.25 and -11.25
        if ~isempty(MatData.class(y).ntr) && length(MatData.class(y).ntr)>1 %% at least 2 error trials
            for i = 1:length(MatData.class(y).ntr)
                ts = MatData.class(y).ntr(i).TS;
                a = find(ts >  MatData.class(y).ntr(i).Sample_onT-time_period & ts < MatData.class(y).ntr(i).Sample_onT);
                erranmrate(i) = length(a)/time_period;
            end
            %     rate(c/2) = ((sum(mrate)+sum(nmrate))/(length(MatData.class(c).ntr)+length(MatData.class(c-1).ntr)))
            errsumrate(y/2) = sum(erranmrate);
            errsumtrialnum(y/2) = length(MatData.class(y).ntr);
            erstdrate(y/2) = std(erranmrate);
            clear erranmrate;
        else
            errsumrate(y/2) = nan;%% if less than 2 error trials,set firing rate at 0
            errsumtrialnum(y/2) = 0;
            erstdrate(y/2) = nan;
            clear erranmrate; 
        end
        err_delay_rate(y/2) = errsumrate(y/2)/errsumtrialnum(y/2);
    end
catch
    erranmrate = nan;
    erstdrate = nan(1,8);
    errsumrate = nan(1,8);
    errsumtrialnum = nan(1,8);
end


%%%%%%%%%%%
errdelay = [err_delay_rate(8) err_delay_rate(7) err_delay_rate(6) err_delay_rate(5) err_delay_rate(1) err_delay_rate(2) err_delay_rate(3) err_delay_rate(4)];
if length(erstdrate) == 7
    errstd = [NaN erstdrate(7) erstdrate(6) erstdrate(5) erstdrate(1) erstdrate(2) erstdrate(3) erstdrate(4)];
elseif length(erstdrate) == 6
    errstd = [NaN NaN erstdrate(6) erstdrate(5) erstdrate(1) erstdrate(2) erstdrate(3) erstdrate(4)];
elseif length(erstdrate) == 5
    errstd = [NaN NaN NaN  erstdrate(5) erstdrate(1) erstdrate(2) erstdrate(3) erstdrate(4)];
elseif length(erstdrate) == 4
    errstd = [NaN NaN NaN NaN erstdrate(1) erstdrate (2) erstdrate (3) erstdrate (4)];
else
    errstd = [erstdrate(8) erstdrate(7) erstdrate(6) erstdrate(5) erstdrate(1) erstdrate(2) erstdrate(3) erstdrate(4)];
end
cordelay = [delay_rate(8) delay_rate(7) delay_rate(6) delay_rate(5) delay_rate(1) delay_rate(2) delay_rate(3) delay_rate(4)];
corstd = [stdrate(8) stdrate(7) stdrate(6) stdrate(5) stdrate(1) stdrate(2) stdrate(3) stdrate(4)];
result = [averagefr errdelay cordelay];
average_err = nanmean(errdelay);
average_err_alltrials = nansum(errsumrate)/sum(errsumtrialnum);

%%to test if a location is higher than 0 degree or lower than 0 degree
%%if it is higher, mark as 1, if lower, marks as -1
exam = zeros(1,8);
for b = 1:8
    if result(b+9)>(1+crit)*result(1);
        exam(b) = 1;
    elseif result(b+9)<(1-crit)*result(1);
        exam(b) = -1;
    else
        exam(b) = 0; %% either equal to the location 17's firing rate, nor in the range decided by criterion
    end
end

frpool = nan(1,36);
stdpool = nan(1,36); %% collect all std number
for c = 1:8
    if exam(c) == 1;
        frpool(c) = cordelay(c);%%extract correct trials firing rate
        frpool(c+8) = errdelay(c);%%extract error trials firing rate
        stdpool(c) = corstd(c); %%extract correct trial's std
        stdpool(c+8) = errstd(c); %%extract error trial's std
    elseif  exam(c) == -1;
        frpool(c+16) = cordelay(c);%%extract correct trials firing rate
        frpool(c+24) = errdelay(c);%%extract error trials firing rate
        stdpool(c+16) = corstd(c); %%extract correct trial's std
        stdpool(c+24) = errstd(c); %%extract error trial's std
    end
end

average_result = nan(1,4);
average_std = nan(1,4);
average_result(1) = nanmean(frpool(1:8));%%FR from correct higher class trials
average_result(2) = nanmean(frpool(9:16));%%FR from error higher class trials
average_result(3) = nanmean(frpool(17:24));%%FR from correct lower class trials
average_result(4) = nanmean(frpool(25:32));%%FR from error lower class trials
average_std(1) = nanmean(stdpool(1:8));
average_std(2) = nanmean(stdpool(9:16));
average_std(3) = nanmean(stdpool(17:24));
average_std(4) = nanmean(stdpool(25:32));

frpool(33:36) = average_result;
stdpool(33:36) = average_std;

frpool_alltrials = nan(2,16);

for d = 1:length(erstdrate)
    frpool_alltrials(1,d) = sumrate(d);
    frpool_alltrials(1,d+8) = sumtrialnum(d);
    frpool_alltrials(2,d) = errsumrate(d);
    frpool_alltrials(2,d+8) = errsumtrialnum(d);
end

frpool_alltrials_reo = nan(2,16);%%reorganize the firing rate database
frpool_alltrials_reo(1,:) = [frpool_alltrials(1,8) frpool_alltrials(1,7) frpool_alltrials(1,6) frpool_alltrials(1,5) frpool_alltrials(1,1) frpool_alltrials(1,2) frpool_alltrials(1,3) frpool_alltrials(1,4) frpool_alltrials(1,16) frpool_alltrials(1,15) frpool_alltrials(1,14) frpool_alltrials(1,13) frpool_alltrials(1,9) frpool_alltrials(1,10) frpool_alltrials(1,11) frpool_alltrials(1,12)];
frpool_alltrials_reo(2,:) = [frpool_alltrials(2,8) frpool_alltrials(2,7) frpool_alltrials(2,6) frpool_alltrials(2,5) frpool_alltrials(2,1) frpool_alltrials(2,2) frpool_alltrials(2,3) frpool_alltrials(2,4) frpool_alltrials(2,16) frpool_alltrials(2,15) frpool_alltrials(2,14) frpool_alltrials(2,13) frpool_alltrials(2,9) frpool_alltrials(2,10) frpool_alltrials(2,11) frpool_alltrials(2,12)];

frpool_alltrials_err = nan(2,32);
for e = 1:8
    if exam(e) == 1;
        frpool_alltrials_err(1,e) = frpool_alltrials_reo(1,e);%%extract correct trials firing rate
        frpool_alltrials_err(2,e) = frpool_alltrials_reo(1,e+8);%%extract correct trials number
        frpool_alltrials_err(1,e+8) = frpool_alltrials_reo(2,e);%%extract error trials firing rate
        frpool_alltrials_err(2,e+8) = frpool_alltrials_reo(2,e+8);%%extract error trials number
        
    elseif  exam(e) == -1;
        frpool_alltrials_err(1,e+16) = frpool_alltrials_reo(1,e);%%extract correct trials firing rate
        frpool_alltrials_err(2,e+16) = frpool_alltrials_reo(1,e+8);%%extract correct trials number
        frpool_alltrials_err(1,e+24) = frpool_alltrials_reo(2,e);%%extract error trials firing rate
        frpool_alltrials_err(2,e+24) = frpool_alltrials_reo(2,e+8);%%extract error trials number
    end
end

fr_alltrials_err = nan(1,4);
fr_alltrials_err(1) = nansum(frpool_alltrials_err(1,1:8))/nansum(frpool_alltrials_err(2,1:8));
fr_alltrials_err(2) = nansum(frpool_alltrials_err(1,9:16))/nansum(frpool_alltrials_err(2,9:16));
fr_alltrials_err(3) = nansum(frpool_alltrials_err(1,17:24))/nansum(frpool_alltrials_err(2,17:24));
fr_alltrials_err(4) = nansum(frpool_alltrials_err(1,25:32))/nansum(frpool_alltrials_err(2,25:32));

excel_output = [result average_result average_std];
% op = [result,233,frpool,233,fr_alltrials_err, exam] % sophia playing
% around
op = [exam];
% % average_result_alltrials = nan(1,4)
% % average_result_alltrials(1) = nansum(
%
% h = figure;
% scatter(location, errdelay, 'r','filled');%%plot error trials
% hold on
% scatter(0, averagefr, 'k','*');
% %plot([-90 90],[average_err average_err],'--');%plot mean error firing rate
% %plot([-90 90],[average_err_alltrials average_err_alltrials],'k');
% plot(location, cordelay, 'b');%%plot correct trials
% title(['last 1s'])