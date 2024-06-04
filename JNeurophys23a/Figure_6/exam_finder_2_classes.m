function [exam, Classes] = exam_finder_2_classes( filename,r1r2)

% load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', filename]);
% time_period = 1;%% test last 0.5/1/3 delay period firing rate
% fr = zeros(1,16);

%% measure the best position
max_class_corr = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', filename]);
if max_class_corr(1) == 1
%     Classes = [3 8];
%     Classes = [13 18];
%     exam=[1 -1];
    if r1r2==1
        Classes=[1 2 3 4 5  6  7  8  9 10]; % for remember 1st
    elseif r1r2==2
        Classes=[11 12 13 14 15 16 17 18 19 20]; % for remember 2nd
%         Classes=[11 16];
    end
    exam   =[1 1 1 1 1 -1 -1 -1 -1 -1];
%     exam   =[1 -1];
elseif max_class_corr(1) ==6
%     Classes = [8 3];
%     Classes = [18 13];
%     exam=[1 -1];
    if r1r2==1
        Classes=[6 7 8 9 10 1  2  3  4  5]; %for remember 1st
    elseif r1r2==2
        Classes=[16 17 18 19 20 11 12 13 14 15]; % for remember 2nd
%         Classes=[16 11];
    end
    exam =  [1 1 1 1 1 -1 -1 -1 -1 -1];
%     exam =  [1 -1];
else
    Classes = [];
    exam=[];
end
% clear MatData


% %%
% 
% %%generate delay rate for different locations from NM trials
% zi = 0;
% for z = Classes
%     try
%         for i = 1:length(MatData.class(z).ntr)
%             ts = MatData.class(z).ntr(i).TS;
%             a = find(ts >  MatData.class(z).ntr(i).Sample_onT-time_period & ts < MatData.class(z).ntr(i).Sample_onT);
%             anmrate(i) = length(a)/time_period;
%         end
%         sumrate(zi+1) = sum(anmrate);
%         sumtrialnum(zi+1) = length(MatData.class(z).ntr);
%     %     stdrate(z/2) = std(anmrate);
%         clear anmrate
%     catch
%         sumrate(zi+1) = nan;
%         sumtrialnum(zi+1) = 0;
% 
%     end
%     zi = zi+1;
% end
% delay_rate = sumrate./sumtrialnum;
% 
% clear ts a MatData
% 
% %%generate delay rate for different locations for error trials
% fn_err = [filename '_err'];
% load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);
% % load([filename '_err']);
% err_delay_rate = nan(1,length(Classes));
% yi=0;
% for y = Classes
%     try
% %     if any(mean([MatData.class(y).ntr.cuerate]))
% %     if ~isempty(MatData.class(y).ntr) && length(MatData.class(y).ntr)>1 %% at least 2 error trials
%         for i = 1:length(MatData.class(y).ntr)
%             ts = MatData.class(y).ntr(i).TS;
%             a = find(ts >  MatData.class(y).ntr(i).Target_onT-time_period & ts < MatData.class(y).ntr(i).Target_onT);
%             erranmrate(i) = length(a)/time_period;
%         end
%         %     rate(c/2) = ((sum(mrate)+sum(nmrate))/(length(MatData.class(c).ntr)+length(MatData.class(c-1).ntr)))
%         errsumrate(yi+1) = sum(erranmrate);
%         errsumtrialnum(yi+1) = length(MatData.class(y).ntr);
% %         erstdrate(y/2) = std(erranmrate);
%         clear erranmrate;
%     catch
%         errsumrate(yi+1) = nan;%% if less than 2 error trials,set firing rate at 0
%         errsumtrialnum(yi+1) = 0;
% %         erstdrate(y/2) = nan;
%         clear erranmrate;
%     end
%     
%     err_delay_rate(yi+1) = errsumrate(yi+1)/errsumtrialnum(yi+1);
%     yi = yi + 1;
% end
% 
% 
% % %%%%%%%%%%%
% errdelay = err_delay_rate;  %[err_delay_rate(8) err_delay_rate(7) err_delay_rate(6) err_delay_rate(5) err_delay_rate(1) err_delay_rate(2) err_delay_rate(3) err_delay_rate(4)];
% % %errstd = [erstdrate(8) erstdrate(7) erstdrate(6) erstdrate(5) erstdrate(1) erstdrate(2) erstdrate(3) erstdrate(4)];
% cordelay = delay_rate; %[delay_rate(8) delay_rate(7) delay_rate(6) delay_rate(5) delay_rate(1) delay_rate(2) delay_rate(3) delay_rate(4)];
% % %corstd = [stdrate(8) stdrate(7) stdrate(6) stdrate(5) stdrate(1) stdrate(2) stdrate(3) stdrate(4)];
% % result = [averagefr errdelay cordelay];
% result = [cordelay(1) errdelay(1) cordelay(2) errdelay(2)];
% % average_err = nanmean(errdelay);
% % average_err_alltrials = nansum(errsumrate)/sum(errsumtrialnum);