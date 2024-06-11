function [MI_matrix_raw, MI_matrix_miu, MI_matrix_delta] = PACmeg_surrogate(cfg,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PACmeg: a function to perform phase amplitude coupling analysis
%
% Author: Robert Seymour (rob.seymour@ucl.ac.uk) April 2020
% Edited: Zhengyang wang Jan 2021
% Requirements:
% - MATLAB 2016b or later
% - Fieldtrip Toolbox
% - PACmeg
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%%
%
% data              = data for PAC (size: trials*time)
% cfg.Fs            = Sampling frequency (in Hz)
% cfg.phase_freqs   = Phase Frequencies in Hz (e.g. [8:1:13])
% cfg.amp_freqs     = Amplitude Frequencies in Hz (e.g. [40:2:100])
% cfg.filt_order    = Filter order used by ft_preproc_bandpassfilter
%
% amp_bandw_method  = Method for calculating bandwidth to filter the
%                   ampltitude signal:
%                        - 'number': +- nHz either side
%                        - 'maxphase': max(phase_freq)
%                        - 'centre_freq': +-2.5*amp_freq
% amp_bandw         = Bandwidth when cfg.amp_bandw_method = 'number';
%
% cfg.method        = Method for PAC Computation:
%                   ('tort','ozkurt','plv','canolty)
%
% cfg.surr_method   = Method to compute surrogates:
%                        - '[]': No surrogates
%                        - 'swap_blocks': cuts each trial's amplitude at
%                        a random point and swaps the order around
%                        - 'swap_trials': permutes phase and amp from
%                        different trials
% cfg.surr_N        = Number of iterations used for surrogate analysis
%
% cfg.mask          = filters ALL data but masks between times [a b]
%                   (e.g. cfg.mask = [100 800]; will
%
% cfg.avg_PAC       = Average PAC over trials ('yes' or 'no')
%
%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%
%
% - MI_matrix_raw   = comodulagram matrix (size: amp*phase)
% - MI_matrix_surr  = surrogate comodulagram matrix (size: surr*amp*phase)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check if Fieldtrip is in the MATLAB path
try
    ft_defaults;
catch
    error('Add Fieldtrip to your MATLAB path');
end

%% Get function inputs
% Get sampling frequency
Fs = ft_getopt(cfg,'Fs',[]);

if isempty(Fs)
    error('Please Specify cfg.Fs');
elseif ~isnumeric(Fs)
    error('cfg.Fs needs to be numeric (e.g. 1000)');
end

% Get phase frequencies
phase_freqs = ft_getopt(cfg,'phase_freqs',[]);

if isempty(phase_freqs)
    error('Please Specify cfg.phase_freqs');
end

% Get amplitude frequencies
amp_freqs = ft_getopt(cfg,'amp_freqs',[max(phase_freqs):2:Fs/2]);

% Get filter order
filt_order = ft_getopt(cfg,'filt_order',4);

% Get amplitude bandwidth method
amp_bandw_method = ft_getopt(cfg,'amp_bandw_method','maxphase');
amp_bandw = ft_getopt(cfg,'amp_bandw',10);

% Get PAC Method
method = ft_getopt(cfg,'method','tort');
% fprintf('Using the %s method for PAC computation\n',method);

% Get Masking
mask = ft_getopt(cfg,'mask',[]);

% Get number of surrogates
n_surr = ft_getopt(cfg,'n_surr', 200);

%% Check inputs

% Check whether the inputs are numbers(!)
if ~floor(phase_freqs) == phase_freqs
    error('Numeric Values ONLY for Phase');
end

if ~floor(amp_freqs) == amp_freqs
    ft_error('Numeric Values ONLY for Amplitude');
end

% Give user a warning if using low-frequencies for phase
if min(phase_freqs) < 7 && filt_order > 3
    ft_warning(['Think about using a lower filter order '...
        '(e.g. cfg.filt_order = 3)']);
end

% If incorrect method abort and warn  the user
if ~any(strcmp({'tort','plv','ozkurt','canolty'},method))
    error(['Specified PAC method ''' method ''' not supported']);
end

% Check whether PAC can be detected
switch amp_bandw_method
    
    case 'number'
        % If the bandwidth is less than the maximum phase frequency...
%         if amp_bandw < max(phase_freqs)
%             
%             error(['You will not be able to detect PAC with this configuration.'...
%                 ' Reduce the phase to ' ...
%                 num2str(amp_bandw) 'Hz, or increase the amplitude bandwidth to '...
%                 num2str(max(phase_freqs)+1) 'Hz']);
%         end
    case 'maxphase'
        % If minimum
        %         if min(amp_freqs) - max(phase_freqs) < max(phase_freqs)
        %             error(['You will not be able to detect PAC with this configuration.'])
        %         end
    case 'centre_freq'
        % If
        if min(amp_freqs)/2.5 < max(phase_freqs)
            try
                low_amp = min(amp_freqs(find(amp_freqs/2.5 > max(phase_freqs))));
            catch
                low_amp = '?';
            end
            
            error(['You will not be able to detect PAC with this configuration.'...
                ' Reduce the phase to ' ...
                num2str(min(amp_freqs)/2.5) 'Hz, or increase the amplitude to '...
                num2str(low_amp) 'Hz']);
        end
end

%% Filter Phase Frequencies & take 'angle'
% disp('Filtering Phase...');


phase_filtered = zeros(size(data,1),length(phase_freqs),length(data));

for phase = 1:length(phase_freqs)
    try
        [filt] = ft_preproc_bandpassfilter(data, Fs,...
            [phase_freqs(phase)-1 phase_freqs(phase)+1],...
            filt_order, 'but');
    catch
        error('Could not filter ... Perhaps try a lower filter order');
    end
    
    
    phase_filtered(:,phase,:) = ft_preproc_hilbert(filt, 'angle');
    
    clear filt
end


%% Filter Amplitude & Take 'abs'
% disp('Filtering Amplitude...');
amp_filtered = zeros(size(data,1),length(amp_freqs),length(data));

for amp = 1:length(amp_freqs)
    
    % Switch based on bandwidth method
    switch amp_bandw_method
        
        case 'number'
            
            if amp == 1
                % fprintf('Bandwidth = %.1fHz\n',amp_bandw);
            end
            
            Af1 = amp_freqs(amp) - amp_bandw;
            Af2 = amp_freqs(amp) + amp_bandw;
            %
        case 'maxphase'
            if amp == 1
                % fprintf('Bandwidth = %.1fHz\n',max(phase_freqs));
            end
            %
            Af1 = amp_freqs(amp) - max(phase_freqs);
            Af2 = amp_freqs(amp) + max(phase_freqs);
            
        case 'centre_freq'
            if amp == 1
                % fprintf('Bandwidth = 2.5* centre amplitude frequency\n')
            end
            
            Af1 = round(amp_freqs(amp) -(amp_freqs(amp)/2.5));
            Af2 = round(amp_freqs(amp) +(amp_freqs(amp)/2.5));
            
            
    end
    
    % Filter
    [filt] = ft_preproc_bandpassfilter(data, Fs,...
        [Af1 Af2],filt_order, 'but');
    
    amp_filtered(:,amp,:) = ft_preproc_hilbert(filt, 'abs');
    
    
    clear filt Af1 Af2
end
%% PAC computation
MI_matrix_raw = zeros(size(mask, 1), length(amp_freqs),length(phase_freqs));
MI_matrix_p = zeros(size(mask, 1), length(amp_freqs),length(phase_freqs));
MI_matrix_surr = zeros(n_surr, size(mask, 1), length(amp_freqs),length(phase_freqs));
for i = 1:n_surr
%     tic
    surr_order = randperm(size(data, 1));
    for phase = 1:length(phase_freqs)
        for amp = 1:length(amp_freqs)
            for i_mask = 1:size(mask, 1)
                %   Keeping phase_data still and permuting amp_data
                phase_data = squeeze(phase_filtered(:, phase,mask(i_mask, 1):mask(i_mask, 2)));
                phase_data = reshape(phase_data, 1, numel(phase_data));
                amp_data = squeeze(amp_filtered(surr_order, amp,mask(i_mask, 1):mask(i_mask, 2)));
                amp_data = reshape(amp_data, 1, numel(amp_data));
                % Switch based on the method of PAC computation
                switch method
                    case 'tort'
                        [MI] = calc_MI_tort(phase_data,amp_data,18);
                        
                    case 'ozkurt'
                        [MI] = calc_MI_ozkurt(phase_data,amp_data);
                        
                    case 'plv'
                        [MI] = cohen_PLV(phase_data,amp_data);
                        
                    case 'canolty'
                        [MI] = calc_MI_canolty(phase_data,amp_data);
                end
                
                % Add to matrix outside the loop
                MI_matrix_surr(i, i_mask, amp, phase) = MI;
            end
        end
    end
%     toc
end

for phase = 1:length(phase_freqs)
    for amp = 1:length(amp_freqs)
        for i_mask = 1:size(mask, 1)
            phase_data = squeeze(phase_filtered(:, phase,mask(i_mask, 1):mask(i_mask, 2)));
            phase_data = reshape(phase_data, 1, numel(phase_data));
            amp_data = squeeze(amp_filtered(:, amp,mask(i_mask, 1):mask(i_mask, 2)));
            amp_data = reshape(amp_data, 1, numel(amp_data));
            switch method
                case 'tort'
                    [MI] = calc_MI_tort(phase_data,amp_data,18);
                    
                case 'ozkurt'
                    [MI] = calc_MI_ozkurt(phase_data,amp_data);
                    
                case 'plv'
                    [MI] = cohen_PLV(phase_data,amp_data);
                    
                case 'canolty'
                    [MI] = calc_MI_canolty(phase_data,amp_data);
            end
            MI_matrix_raw(i_mask, amp, phase) = MI;
            [temp_miu, temp_delta] = normfit(squeeze(MI_matrix_surr(:, i_mask, amp, phase)));
            MI_matrix_miu(i_mask, amp, phase) = temp_miu;
            MI_matrix_delta(i_mask, amp, phase) = temp_delta;
%             MI_matrix_p(i_mask, amp, phase) = 1-normcdf(abs((MI - temp_miu)/temp_delta)); % Two-tail
        end
    end
end
end
