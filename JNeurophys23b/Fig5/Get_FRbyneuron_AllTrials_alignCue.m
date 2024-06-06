function [FR_temp1, ntrs_temp] = Get_FRbyneuron_AllTrials_alignCue(datain,class_num,epoch_start_in,epoch_end_in)
%return the mean firing rate of a certain epoch of each neuron
%for ODR 2012 ver.
%Junda Zhu, 20230619
try
    MatData = datain;
    nTS1 = [];
    m_counter1 = 0;

    epoch_start = epoch_start_in;
    epoch_end = epoch_end_in; % the certain time period according to the cue onset

    if ~isempty(MatData) && class_num <= length(MatData)
        for m1 = 1:length(MatData{class_num})
            if ~isempty([MatData{class_num}.Cue_onT])
                try
                    TS=[];
                    TS = MatData{class_num}(m1).TS-MatData{class_num}(m1).Cue_onT;
                    nTS = length(find(TS>=epoch_start & TS< epoch_end));
                    nTS1 = [nTS nTS1];
                    m_counter1 = m_counter1 + 1;
                catch
                end
            end
        end
        ntrs1 = m_counter1;
    else
        disp('Empty MatData File!!!');
    end

    ntrs_temp = ntrs1;
    FR_temp1 = mean(nTS1/(epoch_end-epoch_start));
catch
    ntrs_temp = nan;
    FR_temp1 = nan;
end