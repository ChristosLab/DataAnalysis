function LFPData = LFP_Data_Extraction_error_edited(folder_name,filename,CHAN)
%
% % Generates MAT file containing LFP data % LFPDATA
% =LFP_DATA('filename',channel) % trials % Edited on 4/25/21 to extract
% error trials only and event timing from diode events - ZW

% persistent trials trials2 lastfn

channel       = CHAN;
trials        = APMReadUserData([folder_name,'/',filename,'.apm']);
diode_channel = 5;
sampling_rate = 500;
delay_thresh = 0.5; %   Monitor delay threshold for diode error 
LFPData       = [];

%% Try to retrieve behavioral data
try
    AllData = [];
    load(fullfile(folder_name, [filename,'.mat']))  %%Loads AllData
catch
    disp(['problem getting behavioral data, exiting ', filename])
    return
end
e_counter = 1;
err_log   = [];
conv_dirs = 1;
trials = trials(1:numel(AllData.trials)); %  Match trial counts between APM and MAT file
%%SCR datafiles contain cue condition
% if length(unique([trials.dir])) > 9  %Convert actual class number (1-36, if present) to cue condition (1-9)
%     conv_dirs = 1;                   %SCR datafiles contain cue condition
% else                                 %ELV datafiles contain class number
%     conv_dirs = 0;
% end

for n = 1:length(trials)
    if conv_dirs
        trials(n).dir = [AllData.trials(n).Class];
    end
end

%%% Local parameters
dirs    = unique([trials.dir]);
numdirs = length(dirs);

%%% Add time points for trial epochs %Datafile contains variable number of stimuli
% if isfield(AllData.ClassStructure, 'frametime')
%     time_frame = 1;
% else  %default time values
%     time_frame = 2;
% end

%% Main Code %%
cl_n = zeros(1, numdirs);
for n = 1:length(trials)
    try
        if (strcmp(AllData.trials(n).Reward, 'No')) && (AllData.trials(n).Statecode > 4) % Saccade error trials only
            %                 if (trials(n).rewarded==1)  %% Process only correct trials
            %         if strcmpi(filename(1:3), 'scr')
            %             diode_loc = trials(n).channels(diode_channel).timestamp;
            %             diode_Hz  = 142.0455;
            %         else
            test_diode = 1;
            tested_channel = [];
            while test_diode
                diode_loc = trials(n).channels(diode_channel).events1;
                diode_Hz  = 25; %% Assuming the sampling rate is 40 kHz!!!!
                if isempty(diode_loc)
                    fprintf('No diode event on channel%d on trial %d in file %s\n', diode_channel, n, filename)
                    test_diode = 1;
                    tested_channel = [tested_channel, diode_channel];
                    rest_of_channels = setxor(1:8, tested_channel);
                    diode_channel = rest_of_channels(1);
                end
                test_diode = 0;
                diode_channel = 5;
            end
            %         end
            %             switch time_frame
            %                 case 1
            %                     % %                       fix_end = AllData.ClassStructure(AllData.trials(n).Class).frametime(1);
            %                     % %                       cue_dur = AllData.ClassStructure(AllData.trials(n).Class).frametime(2);
            %                     % %                       cuedelay_dur = AllData.ClassStructure(AllData.trials(n).Class).frametime(3);
            %                     % %                       sample_dur = AllData.ClassStructure(AllData.trials(n).Class).frametime(4);
            %                     sampledelay_dur    = AllData.ClassStructure(AllData.trials(n).Class).frametime(5);
            %                     if sampledelay_dur == 0
            %                         diode_test     = 2;
            %                     else
            %                         diode_test     = 3;
            %                     end
            %                 case 2
            %                     fix_duration = AllData.parameters.fixationDuration;
            %                     cue_duration = AllData.parameters.stimulusDuration;
            %                     del_duration = AllData.parameters.delayDuration;
            %                     tar_duration = AllData.parameters.targetDuration;
            %                     % %                       cuedelay_dur      = 1.5;
            %                     % %                       sample_dur        = 0.5;
            %                     % %                         sampledelay_dur = 0.5;
            %                     diode_test   = 3;
            %             end
            syncs = diode_loc*diode_Hz*1e-6;  %Actual diodes
            CueOnT = syncs(1);
            beh_CueOnT = AllData.trials(n).CueOn - AllData.trials(n).time; %    CueOnT recorded on 
            if abs(CueOnT - beh_CueOnT) > delay_thresh
                fprintf('Diode error on trial %d in file %s\n', n, filename)
                continue
            end
            %             if length(syncs) >= diode_test    %Expected diodes
            Class_var       = trials(n).dir;
            cl_n(find(dirs == Class_var)) = cl_n(find(dirs == Class_var)) + 1;
            ntr             = cl_n(find(dirs == Class_var));
            LFPData.class(Class_var).ntr(ntr).Cue_onT    = CueOnT;
            %                 LFPData.class(Class_var).ntr(ntr).Cue_onT    = syncs(1);
            %                 switch length(syncs)
            %                     case 4
            %                         LFPData.class(Class_var).ntr(ntr).Cue_onT    = syncs(1);
            %                         LFPData.class(Class_var).ntr(ntr).Sample_onT = syncs(2);
            %                         LFPData.class(Class_var).ntr(ntr).target_onT = syncs(3);
            %                         LFPData.class(Class_var).ntr(ntr).Reward_onT = syncs(4);
            %                     case 3
            %                         LFPData.class(Class_var).ntr(ntr).Cue_onT    = syncs(1);
            %                         LFPData.class(Class_var).ntr(ntr).Sample_onT = syncs(2);
            %                         LFPData.class(Class_var).ntr(ntr).target_onT = syncs(3);
            %                     case 2
            %                         LFPData.class(Class_var).ntr(ntr).Cue_onT    = syncs(1);
            %                         LFPData.class(Class_var).ntr(ntr).target_onT = syncs(2);
            %                     case 1
            %                         LFPData.class(Class_var).ntr(ntr).Cue_onT    = syncs(1);
            %                 end
            %% Fill in LFPData structure
            if isfield(trials(n).channels(diode_channel), 'events1')
                LFPData.class(Class_var).ntr(ntr).Syncs     = length(trials(n).channels(diode_channel).events1);
            else
                LFPData.class(Class_var).ntr(ntr).Syncs     = nan;
            end
            LFPData.class(Class_var).ntr(ntr).Trial_Num = n;
            
            if trials(n).channels(channel).lfp_time_calibration > 1
                LFPData.class(Class_var).ntr(ntr).Sampling_Rate = 1000000/trials(n).channels(channel).lfp_time_calibration;
            else
                LFPData.class(Class_var).ntr(ntr).Sampling_Rate = 500;
            end
%             if (strcmp(AllData.trials(n).Reward,'Yes')) && (trials(n).rewarded==1) %Match = 1, Non-Match = 0;
%                 LFPData.class(Class_var).ntr(ntr).rewarded = 1;  %In some classes match or nonmatch may not appear at all
%             else (strcmp(AllData.trials(n).Reward,'No')) && (trials(n).rewarded==0)
%                 LFPData.class(Class_var).ntr(ntr).rewarded = 0;
%             end
            
            if isfield(trials(n).channels(channel), 'LFP')
                LFPData.class(Class_var).ntr(ntr).LFP = trials(n).channels(channel).LFP;
            end
%             if isfield(trials(n).channels(channel), 'spikes')
%                 LFPData.class(Class_var).ntr(ntr).spikes = trials(n).channels(channel).spikes;
%             end
%             if isfield(trials(n).channels(channel), 'timestamp')
%                 LFPData.class(Class_var).ntr(ntr).timestamp = trials(n).channels(channel).timestamp;
%             end
%         else
%             err_log(e_counter,:).err = ['Unreward trial ', num2str(n)];
%             e_counter                = e_counter + 1;
%             disp(['wrong diode number on trial ' num2str(n)])
        end
        %         else
        %             err_log(e_counter,:).err = ['wrong diode number on trial ', num2str(n)];
        %             e_counter                = e_counter + 1;
        %             disp(['wrong diode number on trial ' num2str(n)])
        %         end
    catch
        err_log(e_counter,:).err = [lasterr, ' on trial ', num2str(n)];
        e_counter                = e_counter + 1;
        disp([lasterr, ' on trial ', num2str(n)])
    end
end
% if ~isempty(err_log)
%     save([filename, '_', num2str(channel), '_log'], 'err_log');
% end
disp('Finished');
end