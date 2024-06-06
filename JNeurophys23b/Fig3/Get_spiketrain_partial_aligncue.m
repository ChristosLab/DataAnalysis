function [spiketrain_temp, ntrs] = Get_spiketrain_partial_aligncue(datain,infoin,range)
% return time stamps in ms. poll all trials
% for ODR task
% 20230607, J Zhu
% datain: neuron data (8 classes);
% infoin: class of choice;
% range: [lo, hi] to set a range to use a subset of the trials/epoches (optional)
class = infoin;
TS_all = {};
m_counter = 0;
for n = 1:length(datain{class})
    try
        m_counter = m_counter + 1;
        TS = [];
        TS = datain{class}(n).TS - datain{class}(n).Cue_onT;
        TS_all{m_counter} = TS*1000;
    catch
    end
end
spiketrain_temp = TS_all';
ntrs = m_counter;
end