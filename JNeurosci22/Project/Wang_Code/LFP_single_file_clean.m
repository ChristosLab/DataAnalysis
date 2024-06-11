function LFP_single_file_clean(in_dir, out_dir, figure_dir)
%%  Path setting
[codlib, ~, ~] = fileparts(mfilename('fullpath'));
chronux_folder = '\External\chronux_2_12\';
ft_folder = '\External\fieldtrip\';
addpath(genpath(fullfile(codlib, chronux_folder)));
addpath(fullfile(codlib, ft_folder));
ft_defaults;
%%  File name parsing options
file_identifier = '*CH*.mat';
excluding_flag = 'catch';
%%  Filtering options. - Balbir Singh
N       = 5; %  2Nth order bandpass filter (See document for BUTTER()); FILTFILT() further results in 4Nth order
flp     = 1;
fhi     = 200;
fs      = 500;
wn      = [flp/(fs/2), fhi/(fs/2)];
[bb,aa] = butter(N,wn,'bandpass');
pl_freq = 60;
notch_w = 1/10;
stop_band = (pl_freq - 1 * (notch_w)):(notch_w):(pl_freq + 1 * (notch_w)); % Powerline noise removal.
%% Spectrogram inspection options - Balbir Singh
params.Fs       = 500; % sampling frequency
params.fpass    = [.5 100]; % band of frequencies to be kept
params.tapers   = [5 7]; % taper parameters
params.pad      = 1; % pad factor for fft
params.err      = [2 0.05];
params.trialave = 0;
movingwin       = [.5 0.05]; %% moving window
%%  Overwite check
assert(isdir(in_dir), 'Input folder does not exit!')
if ~isfolder(out_dir)
    mkdir(out_dir)
elseif ~numel(dir(out_dir)) == 2 %  Check folder empty or not
    prompt = sprintf(...
        'Output directory:\n%s\nis not empty.\nDo you want to overwrite? Y/N: ', ...
        out_dir...
        );
    str = input(prompt,'s');
    ow_flag = 0;
    if (str == 'Y' || str == 'y')
        ow_flag = 1;
    end
end
%%  Get file list
flist_in = dir(fullfile(in_dir, file_identifier));
flist_out = dir(fullfile(out_dir, file_identifier));
[c_fname, i_in] = setdiff({flist_in.name}, {flist_out.name});
% 
ex_list = []; % Files to exlude based on excluding_flag
for i = 1:numel(c_fname)
    if regexp(c_fname{i}, excluding_flag)
        ex_list = [ex_list, i];
    end
end
i_in(ex_list) = [];
flist_in = flist_in(i_in);
%%  File iterations
for i = 1:numel(flist_in)
    data_ = importdata(fullfile(flist_in(i).folder, flist_in(i).name));
    mad_ = get_LFP_sd_dist([data_.class.ntr]);
        for cls = 1:length([data_.class])
            trials = [data_.class(cls).ntr];
            for tr = 1:length(trials)
                tem = trials(tr).LFP';
                filtered_ = filtfilt(bb,aa,tem);
                [filtered_] = ft_preproc_dftfilter(filtered_,fs, stop_band);  %% power line noise removal
                [pow,t,f,~] = mtspecgramc(filtered_, movingwin, params);
                T1 = linspace(1/fs, numel(filtered_)/fs, numel(filtered_));
                T2 = linspace(T1(1), T1(end), numel(t));
                filtered_mad_ = median(abs(filtered_ - median(filtered_)));
                fig_ = figure('Units', 'normalized', 'Position', [.2 .2 .5 .5]);
                subplot(4, 3, [1:2, 4:5]);
                plot(T1, filtered_/filtered_mad_);
                subplot(4, 3, [7:8, 10:11]);
                pcolor(T2,f,pow');shading flat;
                subplot(4, 3, [6, 9]);
                boxplot(mad_)
                hold on
                plot([0, -2], [0, 0] + filtered_mad_, '--k')
                title(filtered_mad_)
                a = '!'; %  Input placeholder
                while not(isnumeric(a)||isletter(a))
                    a = input('Is is trial good (letter GOOD/number BAD)? ','s');
                    a = a(1);
                    while isempty(a)
                        a = input('Is is trial good (letter GOOD/number BAD)? ','s');
                        a = a(1);
                    end
                end
                if isletter(a)
                    LFPDATA.class(cls).ntr(tr).good = 1; %    GOOD trial
                elseif isnumeric(a)
                    LFPDATA.class(cls).ntr(tr).good = 0; %    BAD trial
                    [~, f_part_] = fileparts(flist_in(i).name);
                    fig_name_ = sprintf('%s_%u_%u', f_part_, cls, tr);
                    print(fullfile(figure_dir, fig_name_), '-dpng');
                    disp('Saved file %s ...', fig_name_)
                end
                close(fig_)
            end
        end
end             
%%
end
function out_ = get_LFP_sd_dist(in_struct)
n_ = numel(in_struct);
out_ = zeros(n_, 1);
parfor i = 1:n_
    out_(i) = median(abs(in_struct(i).LFP - median(in_struct(i).LFP)));
end    
end