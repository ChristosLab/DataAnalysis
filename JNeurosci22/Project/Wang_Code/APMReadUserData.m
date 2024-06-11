function apmdata = APMReadUserData(filename )
% apmdata = APMReadUserData( filename )
%  Reads data from an APM data file
%  Input arguments:
%   filename    - input data file
%  Output:
%   apmdata - a structure containing all the data
%          in the file. Type "apmdata" at matlab command prompt to
%          list all the memebers. When using gated recording mode, the
%          function returns a vector of structures, each trial being
%          stored in a different element of the vector.
%%  Saving direction and reward outcome in trials

MSG_NACC_WAVEFORM           = 00;
MSG_ACC_WAVEFORM	        = 01;
MSG_NOTRIAL_NACC_WAVEFORM	= 02;
MSG_NOTRIAL_ACC_WAVEFORM	= 03;
MSG_TIMESTAMP				= 04;
MSG_START_OF_TRIAL			= 05;
MSG_END_OF_TRIAL			= 06;
MSG_EXT_EVENT				= 07;
MSG_TRIGGER_LEVEL			= 08;
MSG_TRIGGER_TIME			= 09;
MSG_TRIGGER_SLOPE			= 10;
MSG_WINDOW_TIME				= 11;
MSG_WINDOW_LOW				= 12;
MSG_WINDOW_HIGH				= 13;
MSG_TIME_CALIBRATION		= 14;
MSG_VOLTAGE_CALIBRATION		= 15;
MSG_SAMPLING_FREQUENCY		= 16;
MSG_REWARD					= 17;
MSG_LFP_TIME_CALIBRATION	= 18;   
MSG_TEMPLATE				= 32;
MSG_CONTINUOUS				= 48;   
MSG_LFP     				= 56;   
MSG_TCPIP_USER0				= hex2dec('10000');
MSG_TCPIP_USER1				= hex2dec('10001');
MSG_TCPIP_USER2				= hex2dec('10002');
MSG_TCPIP_USER3				= hex2dec('10003');
MSG_TCPIP_USER4				= hex2dec('10004');
MSG_TCPIP_USER5				= hex2dec('10005');
MSG_TCPIP_USER6				= hex2dec('10006');
MSG_TCPIP_USER7				= hex2dec('10007');
MSG_TCPIP_USER256           = hex2dec('10100');
MSG_TCPIP_USER257           = hex2dec('10101');

% Switches
increment_trial_on_filename          = 0;
increment_trial_on_first_trial_start = 1;
increment_trial_on_time_reversal     = 0;

if (nargin ~= 1)
    error('The function requires one input argument ...\n Please type "help APMRead" for a description of the function');
end
if isempty(filename)
    [filename, pathname]  = uigetfile('*.apm', 'APM Data File (*.apm)');
    filename              = strcat(pathname, filename);
end

% The unit is coded on bits 8 to 11 when using template matching or on bit 0 when using the acceptance window
%% Initialization
nchan   = 16;                           % Default number of channels
nunits  = 4;                            % Default number of units per channel
trials  = [];
ntr     = zeros(1,nchan);               % trial number; must keep track of it on each channel individually, since data may be interleaved between channels.
nts     = zeros(1,nchan);               % number of timestamps for each channel
last_ts = zeros(1,nchan);               % last timestamp for each channel
nwf     = zeros(1,nchan);               % number of waveforms for each channel
nevt    = zeros(1,nchan);               % number of external events for each channel
nevt1   = zeros(1,nchan);               % number of reward/validate events for each channel
vcal    = repmat(1.41/32767.0,1,nchan); % default voltage calibration
tcal    = repmat(1.0/48000.0,1,nchan);  % default time calibration
lfptcal = repmat(1.0/500.0,1,nchan);    % default time calibration
current_filename = [];                  % A user message sent through TCPIP, telling the filename of the current trial descriptor and data

fdat  = fopen(filename,'rb');
if (fdat~=-1)
    disp(sprintf('Processing %s ...',filename));
    while (~feof(fdat))
        m_code=fread(fdat,1,'uint32');                         % Message type
        if (~feof(fdat))
            m_channel = fread(fdat,1,'uint32');                % Message channel & unit
            m_unit    = bitshift(m_channel,-16);               % Upper word of the channel may contain the unit for certain messages
            m_channel = bitand(m_channel,hex2dec('0000FFFF')); % Lower word of the channel contains the channel
            m_length  = fread(fdat,1,'uint32');                % Message length
            %disp([m_code m_channel m_length]);
            %             Deal here with all messages ...
            switch (m_code)
                case MSG_START_OF_TRIAL
                    m_data = fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    %disp(sprintf('Channel %d, Start-of-Trial at %d !',m_channel,m_data));
                    nts(m_channel)     = 0;
                    nwf(m_channel)     = 0;
                    last_ts(m_channel) = 0;
                    nevt(m_channel)    = 0;      % Reset the event count
                    nevt1(m_channel)   = 0;      % Reset the event count
                    ntr(m_channel)     = ntr(m_channel)+1;
                    % Set current trial parameters that are available so far ...
                    if (~isempty(current_filename))
                        trials(ntr(m_channel)).filename = current_filename;   % The current filename message is sent BEFORE the start-of-trial message
                    end
                    trials(ntr(m_channel)).channels(m_channel).voltage_calibration  = vcal(m_channel);
                    trials(ntr(m_channel)).channels(m_channel).time_calibration     = tcal(m_channel);
                    trials(ntr(m_channel)).channels(m_channel).lfp_time_calibration = lfptcal(m_channel); 
                    for i=1:nunits
                        trials(ntr(m_channel)).channels(m_channel).template(i).unit         = i;
                        try
                            trials(ntr(m_channel)).channels(m_channel).template(i).waveform = templates(m_channel).unit(i).waveform;
                        catch
                            trials(ntr(m_channel)).channels(m_channel).template(i).waveform = [];
                        end
                    end
                case MSG_END_OF_TRIAL
                    m_data = fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    %disp(sprintf('Channel %d, End-of-Trial at %d !',m_channel,m_data));
                case MSG_VOLTAGE_CALIBRATION
                    %disp('Voltage Calibration record ...');
                    vcal(m_channel) = fread(fdat,m_length,'float32');
                case MSG_TIME_CALIBRATION
                    %disp('Time Calibration record ...');
                    tcal(m_channel)=fread(fdat,m_length,'float32');
                case MSG_EXT_EVENT
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    %disp('External Event record ...');
                    if (ntr(m_channel)>0)
                        nevt(m_channel) = nevt(m_channel)+1;
                        trials(ntr(m_channel)).channels(m_channel).events(nevt(m_channel)) = m_data;
                    end
                case MSG_REWARD
                    m_data = fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    %disp('Reward Event record ...');
                    if (ntr(m_channel)>0)
                        nevt1(m_channel) = nevt1(m_channel)+1;
                        trials(ntr(m_channel)).channels(m_channel).events1(nevt1(m_channel)) = m_data;
                    end
                case MSG_LFP_TIME_CALIBRATION 
                    %disp('LFP Time Calibration record ...');
                    lfptcal(m_channel)=fread(fdat,m_length,'float32');
                    try
                        if (ntr(m_channel)>0)
                            trials(ntr(m_channel)).channels(m_channel).lfp_time_calibration=tcal(m_channel);
                        end
                    catch
                        disp(num2str(m_channel));
                    end
                case MSG_TIMESTAMP                     %disp(sprintf(' Channel %d, %d timestamps ...',m_channel,m_length/2));
                    m_data=double(fread(fdat,m_length,'uint32'));    % Read message data, a 32-bit integer
                    if (ntr(m_channel)>0)
                        for i=1:(m_length/2)     % disp(sprintf('  %d : %x',i,m_data(2*(i-1)+1)));
                            nts(m_channel) = nts(m_channel)+1;
                            ts             = (m_data(2*(i-1)+2));
                            trials(ntr(m_channel)).channels(m_channel).timestamp(nts(m_channel)) = ts;
                            if (ts<last_ts(m_channel))
                                disp(sprintf('Time reversal in timestamps on trial in d%',ntr(m_channel)));
                                trials(ntr(m_channel)).channels(m_channel).reversal = nts(m_channel);
                                % Emulate a start of trial. the next  8 lines added on 2007 feb 04 according to steph' script 'APMReadUserData
                                nts(m_channel)   = 1;
                                nwf(m_channel)   = 0;
                                nevt(m_channel)  = 0;      % Reset the event count
                                nevt1(m_channel) = 0;      % Reset the event count
                                ntr(m_channel)   = ntr(m_channel)+1;
                                trials(ntr(m_channel)).channels(m_channel).voltage_calibration  = vcal(m_channel);
                                trials(ntr(m_channel)).channels(m_channel).time_calibration     = tcal(m_channel);
                                trials(ntr(m_channel)).channels(m_channel).lfp_time_calibration = lfptcal(m_channel);
                            end
                            last_ts(m_channel) = ts;
                        end
                    end
                case {MSG_NACC_WAVEFORM,MSG_ACC_WAVEFORM,MSG_NOTRIAL_NACC_WAVEFORM,MSG_NOTRIAL_ACC_WAVEFORM}
                    %disp(sprintf(' Waveform, %d samples ...',m_length-1));
                    % if unit=0 but spike is accepted, that means we are
                    % using the window discriminator instead of template
                    % matching
                    if ((~m_unit) & ((m_code==MSG_ACC_WAVEFORM) | (m_code==MSG_NOTRIAL_ACC_WAVEFORM)))
                        m_unit=1;
                    end
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    if (ntr(m_channel)>0)
                        nwf(m_channel)=nwf(m_channel)+1; % disp(sprintf('%d: %s %d',nwf(m_channel), current_filename, m_length));
                        trials(ntr(m_channel)).channels(m_channel).spikes(nwf(m_channel)).unit      = m_unit;
                        trials(ntr(m_channel)).channels(m_channel).spikes(nwf(m_channel)).timestamp = m_data(1);
                        trials(ntr(m_channel)).channels(m_channel).spikes(nwf(m_channel)).waveform  = m_data(2:m_length);
                    end
                case MSG_CONTINUOUS % save lfp   %disp(sprintf(' Continuous data packet, %d samples ...',m_length-1));
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    if (ntr(m_channel)>0)
                        try
                            trials(ntr(m_channel)).channels(m_channel).continuous = cat(1,trials(ntr(m_channel)).channels(m_channel).continuous,m_data(2:end));
                        catch
                            trials(ntr(m_channel)).channels(m_channel).continuous = m_data(2:m_length);
                        end
                    end
                case MSG_TEMPLATE  %disp(sprintf(' Waveform, %d samples ...',m_length-1));
                    m_data = fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    templates(m_channel).unit(m_unit).waveform = m_data(2:m_length);
                    %  disp(sprintf('%d: %s %d',nwf(m_channel), current_filename, m_length));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).template(m_unit).waveform=templates(m_channel).unit(m_unit).waveform;
                    end
                case MSG_LFP                    %disp(sprintf(' LFP data packet, %d samples ...',m_length-1));
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    %disp(sprintf(' LFP data packet, %d samples at %d samples ...',m_length-1,m_data(1)));
                    if (ntr(m_channel)>0)
                        try
                            trials(ntr(m_channel)).channels(m_channel).LFP = cat(1,trials(ntr(m_channel)).channels(m_channel).LFP,m_data(2:end));
                        catch
                            trials(ntr(m_channel)).channels(m_channel).LFP = m_data(2:m_length);
                        end
                    end
                case MSG_TCPIP_USER1
                    m_data = fread(fdat,m_length*4,'char');    % Read message data, a string in this case
                    current_filename = char(m_data.');
                    disp(sprintf(' Trial Filename : %s ...',current_filename));
                    % Deal with TCPIP user message 256
                case MSG_TCPIP_USER256
                    m_data = fread(fdat,m_length,'uint32');    % Read message data, a string in this case
                    try % disp(sprintf(' Trial Number : %d ...',double(m_data(1))));
                        tNumber = double(m_data(1)); % disp(sprintf('  Trial Type  : %d ...',double(m_data(2))));
                        tType   = double(m_data(2))/100;
                        trials(tNumber).dir = tType;
                        %                         disp(sprintf('  Target Dir  : %.2f ...',double(m_data(3))/100));
                        %                         disp(sprintf('  Target Ecc  : %.2f ...',double(m_data(4))/100));
                    catch
                        lasterr;
                    end
                case MSG_TCPIP_USER257
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a string in this case
                    try % disp(sprintf(' Trial Number : %d ...',double(m_data(1))));
                        tNumber = double(m_data(1)); % disp(sprintf('  Rewarded    : %d ...',double(m_data(2))));
                        tReward = double(m_data(2));
                        trials(tNumber).rewarded = tReward;
                    catch
                        lasterr;
                    end
                otherwise
                    %disp(sprintf('Unknown message code %d on channel %d, length %d !',m_code,m_channel,m_length));
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
            end
        end
    end
    fclose(fdat);
else
    disp(sprintf('Could not open %s ...',filename));
end
apmdata=trials;