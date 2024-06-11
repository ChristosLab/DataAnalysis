function Neuron_Data_PSTH_PROSaccade_BestCueFrom8Locs(Excel_neurons)
% for prosac task
% find best cue location for ProSac
% 28-Feb-2012, xzhou

[Neurons_num Neurons_txt] = xlsread([Excel_neurons '.xlsx']);
warning off MATLAB:divideByZero
Neurons = [Neurons_txt(:,1) num2cell(Neurons_num(:,1))];

% Best_Cue = Get_Maxes(Neurons);
opp_index = [5 6 7 8 1 2 3 4 9];
% for n = 1:length(Best_Cue)
%     Opp_Cue(n) = opp_index(Best_Cue(n));
% end

RealBest_Cue = Get_RealMaxes(Neurons);
for n = 1:length(RealBest_Cue)
    RealOpp_Cue(n) = opp_index(RealBest_Cue(n));
end

for n = 1:length(Neurons)
    
    Profilename = [Neurons{n,1},'_',num2str(Neurons{n,2})];   
%     try
%         Propsth1(n,:) = Get_PsthM(Profilename,Best_Cue(n));
%     catch
%     end
%     try
%         Propsth2(n,:) = Get_PsthM(Profilename,Opp_Cue(n));
%     catch
%     end
    try
        Propsth3(n,:) = Get_PsthM(Profilename,RealBest_Cue(n));
    catch
    end
    try
        Propsth4(n,:) = Get_PsthM(Profilename,RealOpp_Cue(n));
    catch
    end
    %         meanvarfix1(n,:)= Neuron_Data_baseline(Antifilename,Best_Cue(n));
    
end


definepsthmax=40;
figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% colors = 'rgb';
% bin_width = 0.05;  % 50 milliseconds bin
% bin_edges_Pro=-1:bin_width:3;
% bins_Pro = bin_edges_Pro+0.5*bin_width;
% 
% hold on
% Propsth1 = mean(Propsth1);
% Propsth2 = mean(Propsth2);
% plot(bins_Pro,Propsth1,'b','LineWidth',1.5);
% plot(bins_Pro,Propsth2,'r','LineWidth',1.5);
% maxpsth3=max([max(Propsth1) max(Propsth2)]);
% line([0 0], [0 60],'color','k')
% line([.5 .5], [0 60],'color','k')
% axis([-1 3 0 definepsthmax+0.2])
% xlabel('Time s')
% ylabel('Firing Rate spikes/s')


colors = 'rgb';
bin_width = 0.05;  % 50 milliseconds bin
bin_edges_Pro=-1:bin_width:3;
bins_Pro = bin_edges_Pro+0.5*bin_width;

hold on
Propsth3 = mean(Propsth3);
Propsth4 = mean(Propsth4);
plot(bins_Pro,Propsth3,'b','LineWidth',1.5);
plot(bins_Pro,Propsth4,'r','LineWidth',1.5);
maxpsth4=max([max(Propsth3) max(Propsth4)]);
line([0 0], [0 60],'color','k')
line([.5 .5], [0 60],'color','k')
line([2 2], [0 60],'color','k')
axis([-1 3 0 definepsthmax+0.2])
xlabel('Time s')
ylabel('Firing Rate spikes/s')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function max_results = Get_Maxes(Neurons)
% max_result(1:length(Neurons),1:3) = NaN;
% for n = 1:length(Neurons)
%     Profilename = [Neurons{n,1},'_',num2str(Neurons{n,2})];
%     Antifilename = [Neurons{n,1}([1:6]),'_2_',num2str(Neurons{n,2})];
%     temp = Neuron_Data_Maxcuerate_ProFrom4LOC(Profilename,Antifilename);
%     max_results(n,1:length(temp)) = temp(1);
% end

function Realmax_results = Get_RealMaxes(Neurons)
Realmax_result(1:length(Neurons),1:3) = NaN;
for n = 1:length(Neurons)
    filename = [Neurons{n,1},'_',num2str(Neurons{n,2})];    
    temp = Neuron_Data_Maxcuerate_ProFrom8LOC(filename);
    Realmax_results(n,1:length(temp)) = temp(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psth_temp = Get_PsthM(filename,class_num)     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%change
load(filename)
bin_width = 0.05;  % 50 milliseconds bin
bin_edges=-.8:bin_width:1.5;
bins = bin_edges+0.5*bin_width;
bin_edges_Pro=-1:bin_width:3;
bins_Pro = bin_edges_Pro+0.5*bin_width;

allTS = [];
m_counter = 0;
if length(MatData.class)>8 %24  % AntiSac
    if ~isempty(MatData)
        for m = 1:length(MatData.class(class_num).ntr)
            if ~isempty(MatData.class(class_num).ntr(m).alignSac_onT)
                try
%                     TS = MatData.class(class_num).ntr(m).TS-MatData.class(class_num).ntr(m).alignR_onT;
%                     TS = MatData.class(class_num).ntr(m).TS-MatData.class(class_num).ntr(m).alignSac_onT;
                    TS = MatData.class(class_num).ntr(m).TS-MatData.class(class_num).ntr(m).Cue_onT;
                    allTS = [allTS TS];
                    m_counter = m_counter + 1;
                catch
                end
            end
        end
        ntrs = m_counter;
    else
        disp('Empty MatData File!!!');
    end
    psth_temp =histc(allTS,bin_edges)/(bin_width*ntrs);
end

if length(MatData.class)== 8  % ProSac
    if ~isempty(MatData)
        for m = 1:length(MatData.class(class_num).ntr)
            try
                TS = MatData.class(class_num).ntr(m).TS-MatData.class(class_num).ntr(m).Cue_onT;
%                 TS = MatData.class(class_num).ntr(m).TS-MatData.class(class_num).ntr(m).alignSac_onT;
                allTS = [allTS TS];
                m_counter = m_counter + 1;
            catch
            end
        end
        ntrs = m_counter;
    else
        disp('Empty MatData File!!!');
    end
    psth_temp =histc(allTS,bin_edges_Pro)/(bin_width*ntrs);
end





