% Make dataset of Depth/layer of PFC neurons in ODR task
% J Zhu, 20230619
%% read target neuron numbers from excel
clearvars
excel_name = 'AllNeurons2016FirstCohort';
sheet_name = 'all';
neuron_wanted = readtable([excel_name '.xlsx'],'sheet',sheet_name);
%% find info and read neuron data from master excel
master_excel_name = 'depth_master';
master_sheet_name = 'all';
neuron_info = readtable([master_excel_name '.xlsx'],'sheet',master_sheet_name);
%% join data
join_t = innerjoin(neuron_wanted,neuron_info);
%% load neural data
neuron_filename = [join_t.Filename, num2cell(join_t.Neuron)];
Best_Cue = Get_Maxes(neuron_filename);
% Best_Cue = Get_Maxes_evoked(neuron_filename);
Best_Cue = Best_Cue(:,1);
Best_Del = Get_Maxes_Del(neuron_filename);
% Best_Del = Get_Maxes_Del_evoked(neuron_filename);
Best_Del = Best_Del(:,1);
Best_Sac = Get_Maxes_Sac(neuron_filename);
Best_Sac = Best_Sac(:,1);
odr_data = cell(length(neuron_filename),8);
for n = 1:length(neuron_filename)
    profilename = [neuron_filename{n,1}(1:6),'_1_',num2str(neuron_filename{n,2})];
    try
        odr_data_temp = get_odr_data(profilename);
        odr_data(n,:) = odr_data_temp;
    catch
        odr_data(n,:) = [];
    end
end
%% save dataset
script = char(fread(fopen([mfilename, '.m'])))';
neuron_info = join_t;
neuron_info.cue = Best_Cue;
neuron_info.del = Best_Del;
neuron_info.sac = Best_Sac;
save("odr_data_depth_20230906_raw_max.mat","odr_data","neuron_info","script");
disp('Finish running')
clearvars
%%%%%%%%
function odr_data_temp = get_odr_data(profilename)
load(profilename)
odr_data_temp = cell(1,8);
for ic = 1:8
    try
        odr_data_temp{ic} = [MatData.class(ic).ntr];
    catch
    end
end
end
%%%%%%%%
function max_results = Get_Maxes_evoked(Neurons)
max_results(1:length(Neurons),1:2) = NaN;
for n = 1:length(Neurons)
    Profilename = [Neurons{n,1}(1:6),'_1_',num2str(Neurons{n,2})];
    try
        temp = Neuron_Data_Maxcuerate_evoked_ProFrom8LOC(Profilename);
        max_results(n,1:length(temp)) = temp(1);
    catch
    end
end
end
%%%%%%%%
function max_results = Get_Maxes_Del_evoked(Neurons)
max_results(1:length(Neurons),1:2) = NaN;
for n = 1:length(Neurons)
    Profilename = [Neurons{n,1}(1:6),'_1_',num2str(Neurons{n,2})];
    %     Antifilename = [Neurons{n,1}([1:6]),'_2_',num2str(Neurons{n,2})];
    try
    temp = Neuron_Data_Maxdelrate_evoked_ProFrom8LOC(Profilename);
    max_results(n,1:length(temp)) = temp(1);
    catch
    end
end
end
function max_results = Get_Maxes(Neurons)
max_results(1:length(Neurons),1:2) = NaN;
for n = 1:length(Neurons)
    Profilename = [Neurons{n,1}(1:6),'_1_',num2str(Neurons{n,2})];
    try
        temp = Neuron_Data_Maxcuerate_ProFrom8LOC(Profilename);
        max_results(n,1:length(temp)) = temp(1);
    catch
    end
end
end
%%%%%%%%
function max_results = Get_Maxes_Del(Neurons)
max_results(1:length(Neurons),1:2) = NaN;
for n = 1:length(Neurons)
    Profilename = [Neurons{n,1}(1:6),'_1_',num2str(Neurons{n,2})];
    %     Antifilename = [Neurons{n,1}([1:6]),'_2_',num2str(Neurons{n,2})];
    try
    temp = Neuron_Data_Maxdelrate_ProFrom8LOC(Profilename);
    max_results(n,1:length(temp)) = temp(1);
    catch
    end
end
end
%%%%%%%%
function max_results = Get_Maxes_Sac(Neurons)
max_results(1:length(Neurons),1:3) = NaN;
for n = 1:length(Neurons)
    Profilename = [Neurons{n,1}(1:6),'_1_',num2str(Neurons{n,2})];
    %     Antifilename = [Neurons{n,1}([1:6]),'_2_',num2str(Neurons{n,2})];
    try
    temp = Neuron_Data_Maxsacrate_ODR_8LOC(Profilename);
    max_results(n,1:length(temp)) = temp(1);
    catch
    end
end
end