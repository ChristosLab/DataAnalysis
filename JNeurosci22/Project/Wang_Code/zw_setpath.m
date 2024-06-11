%%  Path setting
project_dir = fileparts(fileparts(mfilename('fullpath')));
code_lib = 'Wang_Code';
fig_lib = 'Wang_Figure';
% fig_lib = fullfile(project_dir,'//Wang_Figure//');
lfp_extract_scripts = 'LFP_extraction';
output_database = 'Wang_Database';
% output_database = fullfile(project_dir,'//Wang_Database//');
antisac_database = fullfile(project_dir, '//LFPDatabychannel_Antisac//'); % Repository to antisaccade data.
prosac_database = fullfile(project_dir, '//LFPDatabychannel_Prosac//'); % Repository to prosaccade data.
neuron_database = fullfile(project_dir, '//Neuron_Data//');
apm_database = fullfile(project_dir, '//APM_Data//');
mat_database = fullfile(project_dir, '//MAT_Data//');
% apm_database = fullfile(project_dir, '//APM_Data//MissingBatch//');
% apm_database_2 = fullfile(project_dir, '//APM_Data//SecondBatch//');
lfp_database       = fullfile(project_dir, '//LFP_Data//');
lfp_database_error = fullfile(project_dir, '//LFP_Data_error//');
lfp_database2      = fullfile(project_dir, '//LFP_Data_temp//');
apm_file_identifier = '*/*.apm';
file_identifier = '*CH*.mat';
apm_files = dir(fullfile(apm_database, apm_file_identifier));
% apm_files_2 = dir(fullfile(apm_database_2, apm_file_identifier));
lfp_files = dir(fullfile(lfp_database, file_identifier));
% antisac_files = dir(fullfile(antisac_database, file_identifier));
% prosac_files = dir(fullfile(prosac_database, file_identifier));
% all_files = [antisac_files; prosac_files];
addpath(fullfile(project_dir, code_lib));
% addpath(genpath(fullfile(project_dir, lfp_extract_scripts)));