%%  LFP extraction Demo
%   This script extracts LFP signal from a single channel from 1 APM file,
%   1 behavior MAT file and 1 neuron file from a matching sessions.
%   Dependencies: 1) LFP_Data_Extraction.m; 
%                 2) match_timestamp.m
%                 3) APMReadUserData.m
% 
% 
%   Customize your file names below:
% 
apm_file = './kenXXX_1.apm';
neuron_file = ''; % Find matching neuron file
% 
% 
[folder, fname, ~] = fileparts(apm_file);
fname = lower(fname);
%   Defines which channels to process. 
%   One possibility is to only process channels with neurons
channels = 1; 
%   LFPData created here uses timestamp information from the behavioral
%   file under the same directory
%   Look into commented lines of LFP_DATA_EXTRACTION to see examples of
%   reading photodiode information in the APM file
load(neuron_file);
LFPData = LFP_Data_Extraction(folder, fname, channels);
%   MATCH_TIMESTAMP extracts timestamp information from a matching neuron
%   file
[LFPData, err_code] = match_timestamp(LFPData, MatData);
save(fullfile(folder, [fname, '_CH', num2str(channels), '.mat']), 'LFPData');
