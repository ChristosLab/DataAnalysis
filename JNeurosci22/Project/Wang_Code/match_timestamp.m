%%
function [LFPData, err_code] = match_timestamp(LFPData, neuron_data)
err_code = 0;
delay_threshold = 0.4;
for i = 1:numel(LFPData.class)
    if isempty(LFPData.class(i).ntr)
        continue
    else
%         [lb, li] = sort([LFPData.class(i).ntr.Trial_Num]);
%         [nb, ni] = sort([neuron_data.class(i).ntr.TrialNumber]);
%         if ~all(lb == nb) % Check if trial numbers match
%             LFPData = [];
%             err_code = 1;
%             return
%         end
%         for j = 1:numel(li)
%             %   Check wheter the timestamp in the neuron file is too
%             %   different from the one from behavior file
%             %   li(j): LFP trial index
%             %   ni(j): neuron trial index
%             if abs(LFPData.class(i).ntr(li(j)).Cue_onT - neuron_data.class(i).ntr(ni(j)).Cue_onT) < delay_threshold
%             LFPData.class(i).ntr(li(j)).Cue_onT = neuron_data.class(i).ntr(ni(j)).Cue_onT;
%             else
%                 LFPData = [];
%                 err_code = 2;
%             end
%         end
        %   Class mismatch
        try
            [b, ni, li] = intersect([neuron_data.class(i).ntr.TrialNumber], [LFPData.class(i).ntr.Trial_Num]);
        catch
            err_code = 1;
            return
        end
        %   No trials match
        if isempty(b)
            err_code = 1;
            return
        %   Some trials do not match
        elseif numel(b) < max(numel([neuron_data.class(i).ntr.TrialNumber]), numel([LFPData.class(i).ntr.Trial_Num]))
            err_code = 3;
            neuron_data.class(i).ntr = neuron_data.class(i).ntr(ni);
            LFPData.class(i).ntr = LFPData.class(i).ntr(li); %  Remove unmatched trials
            [b, ni, li] = intersect([neuron_data.class(i).ntr.TrialNumber], [LFPData.class(i).ntr.Trial_Num]);
        end
        for j = 1:numel(b)
            %   Check wheter the timestamp in the neuron file is too
            %   different from the one from behavior file
            %   li(j): LFP trial index
            %   ni(j): neuron trial index
%             i
%             li(j)
%             ni(j)
%             if abs(LFPData.class(i).ntr(li(j)).Cue_onT - neuron_data.class(i).ntr(ni(j)).Cue_onT) < delay_threshold
            LFPData.class(i).ntr(li(j)).Cue_onT = neuron_data.class(i).ntr(ni(j)).Cue_onT;
            if isfield(neuron_data.class(i).ntr(ni(j)), 'Saccade_onT')
                LFPData.class(i).ntr(li(j)).Saccade_onT = neuron_data.class(i).ntr(ni(j)).Saccade_onT;
            end
%                 err_code = 2;
%             end
        end

    end
end
end