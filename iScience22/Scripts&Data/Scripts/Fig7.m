%% This script for Figure 3A_B
clc;  clear all; close all;
default_dir = pwd;
% 
Task     = {'Spatial','Feature'};
Training = {'PRETRAINING','POSTTRAINING'};
Area     = {'Mid-Dorsal','Posterior-Dorsal','Posterior-Ventral'};


Time = importdata(['ForLinePlots/TimeWhole.mat']);
T    = linspace(-2,5,length(Time));
T    = T(8:end);
%% Posterior Dorsal
AfterTraining_Active_Alpha  = importdata(['ForLinePlots/AfterTraining_ActiveSpatial_Alpha_InducedPow_ITI.mat']);
AfterTraining_Passive_Alpha = importdata(['ForLinePlots/AfterTraining_PassiveSpatial_Alpha_InducedPow_ITI.mat']);

AfterTraining_Active_Beta  = importdata(['ForLinePlots/AfterTraining_ActiveSpatial_Beta_InducedPow_ITI.mat']);
AfterTraining_Passive_Beta = importdata(['ForLinePlots/AfterTraining_PassiveSpatial_Beta_InducedPow_ITI.mat']);

AfterTraining_Active_GammaLo  = importdata(['ForLinePlots/AfterTraining_ActiveSpatial_GammaLo_InducedPow_ITI.mat']);
AfterTraining_Passive_GammaLo = importdata(['ForLinePlots/AfterTraining_PassiveSpatial_GammaLo_InducedPow_ITI.mat']);

AfterTraining_Active_GammaHi = importdata(['ForLinePlots/AfterTraining_ActiveSpatial_GammaHi_InducedPow_ITI.mat']);
AfterTraining_Passive_GammaHi= importdata(['ForLinePlots/AfterTraining_PassiveSpatial_GammaHi_InducedPow_ITI.mat']);


fig = figure;
set(fig, 'Position', [30 40 850 650]);
subplot(4,1,1); hold on;
ForLineFigure(AfterTraining_Active_Alpha(:,8:end),T,'b'); 
ForLineFigure(AfterTraining_Passive_Alpha(:,8:end),T,'r');
xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;%ylim([-15 15]); 
 title('PD-8-14 (Hz)');hold off;
 

subplot(4,1,2); hold on;
ForLineFigure(AfterTraining_Active_Beta(:,8:end),T,'b'); 
ForLineFigure(AfterTraining_Passive_Beta(:,8:end),T,'r');
xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;%ylim([-15 15]); 
title('PD-16-32 (Hz)');hold off;


subplot(4,1,3); hold on
ForLineFigure(AfterTraining_Active_GammaLo(:,8:end),T,'b'); 
ForLineFigure(AfterTraining_Passive_GammaLo(:,8:end),T,'r');
xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; %ylim([-6 1]);
title('PD-33-64 (Hz)');hold off;


subplot(4,1,4); hold on
ForLineFigure(AfterTraining_Active_GammaHi(:,8:end),T,'b'); 
ForLineFigure(AfterTraining_Passive_GammaHi(:,8:end),T,'r');
xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)');axis tight; %ylim([-6 1]);
title('PD-65-100 (Hz)');hold off;

set(gcf,'renderer','Painters');set(gcf,'name',['AfterTrainingActivePassiveSpatialTask_Induced']);
% print(fig,['Fig7.eps'],'-depsc');
% saveas(fig,['Fig7.fig'],'fig' );
% close all;  
