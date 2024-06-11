%% This script for Figure Error trials 
clc;  clear all; close all;
default_dir = pwd;
% 
Task     = {'Spatial','Feature'};
Training = {'PRETRAINING','POSTTRAINING'};
Area     = {'Mid-Dorsal','Posterior-Dorsal','Posterior-Ventral'};


Time = importdata(['TimeWhole.mat']);
T    = linspace(-2,5,length(Time));
T    = T(8:end);
%% Spatial task Posterior Dorsal
ERR_SPD_Alpha   = importdata(['ForLinePlots/',Task{1},'_',Area{2},'_Alpha_Error_InducedPow_ITI.mat']);
ERR_SPD_GammaLo = importdata(['ForLinePlots/',Task{1},'_',Area{2},'_GammaLo_Error_InducedPow_ITI.mat']);
                           
Pos_SPD_Alpha   = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{2},'_Alpha_InducedPow_ITI.mat']);
Pos_SPD_GammaLo = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{2},'_GammaLo_InducedPow_ITI.mat']);


fig = figure;
set(fig, 'Position', [30 40 850 650]);

subplot(2,3,1); hold on;
ForLineFigure(ERR_SPD_Alpha(:,8:end),T,'b'); 
ForLineFigure(Pos_SPD_Alpha(:,8:end),T,'m');
xline(-1,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; ylim([-5.6 10])
 title('PD-8-14 (Hz)');hold off; 
clear ERR_SPD_Alpha Pos_SPD_Alpha;

subplot(2,3,4); hold on;
ForLineFigure(ERR_SPD_GammaLo(:,8:end),T,'b'); 
ForLineFigure(Pos_SPD_GammaLo(:,8:end),T,'m');
xline(-1,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; 
 title('PD-33-64 (Hz)');hold off;
clear ERR_SPD_GammaLo Pos_SPD_GammaLo;





% %% Mid-Dorsal
ERR_SMD_Alpha   = importdata(['ForLinePlots/',Task{1},'_',Area{1},'_Alpha_Error_InducedPow_ITI.mat']);
ERR_SMD_GammaLo = importdata(['ForLinePlots/',Task{1},'_',Area{1},'_GammaLo_Error_InducedPow_ITI.mat']);

                           
Pos_SMD_Alpha   = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{1},'_Alpha_InducedPow_ITI.mat']);
Pos_SMD_GammaLo = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{1},'_GammaLo_InducedPow_ITI.mat']);
                          
subplot(2,3,2); hold on;
ForLineFigure(ERR_SMD_Alpha(:,8:end),T,'b'); 
ForLineFigure(Pos_SMD_Alpha(:,8:end),T,'m');
xline(-1,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;
 title('MD-8-14 (Hz)');hold off;
clear ERR_SMD_Alpha Pos_SMD_Alpha;



subplot(2,3,5); hold on;
ForLineFigure(ERR_SMD_GammaLo(:,8:end),T,'b'); 
ForLineFigure(Pos_SMD_GammaLo(:,8:end),T,'m');
xline(-1,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;
 title('MD-33-64 (Hz)');hold off;
clear ERR_SMD_GammaLo Pos_SMD_GammaLo;




%% Posterior-Ventral
ERR_SPV_Alpha   = importdata(['ForLinePlots/',Task{1},'_',Area{3},'_Alpha_Error_InducedPow_ITI.mat']);
ERR_SPV_GammaLo = importdata(['ForLinePlots/',Task{1},'_',Area{3},'_GammaLo_Error_InducedPow_ITI.mat']);

                           
Pos_SPV_Alpha   = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{3},'_Alpha_InducedPow_ITI.mat']);
Pos_SPV_GammaLo = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{3},'_GammaLo_InducedPow_ITI.mat']);

                           
                           
subplot(2,3,3); hold on;
ForLineFigure(ERR_SPV_Alpha(:,8:end),T,'b'); 
ForLineFigure(Pos_SPV_Alpha(:,8:end),T,'m');

xline(-1,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;
 title('PV-8-14 (Hz)');hold off;
clear ERR_SPV_Alpha Pos_SPV_Alpha;



subplot(2,3,6); hold on;
ForLineFigure(ERR_SPV_GammaLo(:,8:end),T,'b'); 
ForLineFigure(Pos_SPV_GammaLo(:,8:end),T,'m');

xline(-1,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;
 title('PV-33-64 (Hz)');hold off; ylim([-2.1 2])
clear ERR_SPV_GammaLo Pos_SPV_GammaLo;


set(gcf,'renderer','Painters');set(gcf,'name',['SpatialTaskError_Induced']);
% print(fig,['Fig6.eps'],'-depsc');
% saveas(fig,['Fig6.fig'],'fig' );
% close all;  
