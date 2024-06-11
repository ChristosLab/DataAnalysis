%%  CWT parameters
fs = 500;
cue_dur = [-1 + 1/fs, 6]; % Signal length 
sac_dur = [-6 + 1/fs, 1];
normalizer = 4;
kernel_flag = 3;
down_sample = 10;
f_range = 2:2:128;
avg_method = 1; %   0 - arithmetic; 1 - geometric
%%
load('FIO069_1_CH124.mat');
%%
LFPData_mod = LFPData;
for i = 1:numel(LFPData.class)
    for j = 1:numel(LFPData.class(i).ntr)
        LFPData_mod.class(i).ntr(j).Cue_onT = LFPData_mod.class(i).ntr(j).Cue_onT + 1;
    end
end
%%
fiona_tfr = zw_repo_cwt_temp_1_30_22(...
    LFPData_mod, ...
    cue_dur, sac_dur, normalizer, ...
    kernel_flag, down_sample, f_range, fs...
    );
