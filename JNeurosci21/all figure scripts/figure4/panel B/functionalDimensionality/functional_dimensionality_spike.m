function [subject_IDs_output, test_runs, winning_models, test_correlations] = functional_dimensionality_spike(wholebrain_all, mask, subject_IDs, full, varargin)
%% ???Copyright 2018, Christiane Ahlheim???
%% This program is free software: you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation, either version 3 of the License, or
%% (at your option) any later version.
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% You should have received a copy of the GNU General Public License
%% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Load betas into wholebrain_all; 
% Set 'sphere' to perform searchlight, else ROI.
% Set 'spmfile' to perform noise normalisation/prewhitening if desired
% Residuals and beta estimates must not include any nans; mask with mask_mat

% output:
% bestn: for each run, the winning dimensionality
% r_outer: for each run, the correlation of the winnning dimensionality-reconstruction with the test data set
% r_alter: for each run, correlation achieved by a full-dimensional reconstruction with the test data set

preproc = inputParser;

addOptional(preproc,'sphere',-1,@isnumeric)
addParameter(preproc,'spmfile','',@ischar);
parse(preproc,varargin{:})


prewhiten = false;

spm_path = preproc.Results.spmfile;
if ~(spm_path == '')
    try
        SPM_file = load(spm_path);
        SPM = SPM_file.SPM;
        prewhiten = true;
    catch
        warning('Cant open the SPM file');
    end
end
    
mask_mat = mask;

n_voxel   = sum(mask_mat(:));
n_subject = length(wholebrain_all);
test_correlations=[];
test_runs=[];
winning_models=[];
subject_IDs_output=[];

if (length(subject_IDs) ~= n_subject)
    error('"subject_IDs" must have the same length as "wholebrain_all"');
end

for i_subject = 1:n_subject

    data = wholebrain_all{i_subject}; 
    data = data(logical(mask_mat(:)),:,:);
    n_sessions = size(data,3);
    
    if prewhiten
        % Load residuals.
        VRes = get_residuals_mat(SPM, NaN, mask_vol);
        res = VRes(:, logical(mask_mat(:)))';
        res = reshape(res, size(res, 1), size(res,2)/n_sessions, n_sessions);
        % data and res should have the same number of voxels and not include any nans.
    else
        res = false;
    end
    
    [subject_ID, test_run, winning_model, test_correlation] = roi_estimate_dim(data, res, subject_IDs(i_subject), full);

    
    test_runs = [test_runs; reshape(test_run, [numel(test_run),1])];
    winning_models = [winning_models; reshape(winning_model, [numel(winning_model),1])];
    test_correlations = [test_correlations; reshape(test_correlation, [numel(test_correlation),1])];
    subject_IDs_output = [subject_IDs_output; repmat(subject_ID, [numel(test_run),1])];
        
end




end
