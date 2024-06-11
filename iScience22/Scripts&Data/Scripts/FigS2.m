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
Pre_SPD_Alpha   = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{2},'_Alpha_InducedPow_ITI.mat']);
Pos_SPD_Alpha   = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{2},'_Alpha_InducedPow_ITI.mat']);

Pre_SPD_Beta    = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{2},'_Beta_InducedPow_ITI.mat']);
Pos_SPD_Beta    = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{2},'_Beta_InducedPow_ITI.mat']);

Pre_SPD_GammaLo = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{2},'_GammaLo_InducedPow_ITI.mat']);
Pos_SPD_GammaLo = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{2},'_GammaLo_InducedPow_ITI.mat']);

Pre_SPD_GammaHi = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{2},'_GammaHi_InducedPow_ITI.mat']);
Pos_SPD_GammaHi = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{2},'_GammaHi_InducedPow_ITI.mat']);


fig = figure;
set(fig, 'Position', [30 40 850 650]);


subplot(4,3,1); hold on;
ForLineFigure(Pre_SPD_Alpha(:,8:end),T,'c'); 
ForLineFigure(Pos_SPD_Alpha(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; ylim([-1 16]); 
 title('PD-8-14 (Hz)');hold off;
clear Pre_SPD_Alpha Pos_SPD_Alpha;

subplot(4,3,4); hold on;
ForLineFigure(Pre_SPD_Beta(:,8:end),T,'c'); 
ForLineFigure(Pos_SPD_Beta(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; ylim([-1 10]); 
title('PD-16-32 (Hz)');hold off;
clear Pre_SPD_Beta Pos_SPD_Beta;

subplot(4,3,7); hold on
ForLineFigure(Pre_SPD_GammaLo(:,8:end),T,'c'); 
ForLineFigure(Pos_SPD_GammaLo(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; ylim([-2.5 2.5]);
title('PD-33-64 (Hz)');hold off;
clear Pre_SPD_GammaLo Pos_SPD_GammaLo;

subplot(4,3,10); hold on
ForLineFigure(Pre_SPD_GammaHi(:,8:end),T,'c'); 
ForLineFigure(Pos_SPD_GammaHi(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)');axis tight; ylim([-1 2.2]);
title('PD-65-100 (Hz)');hold off;
clear Pre_SPD_GammaHi Pos_SPD_GammaHi;

%% Mid-Dorsal
Pre_SMD_Alpha   = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{1},'_Alpha_InducedPow_ITI.mat']);
Pos_SMD_Alpha   = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{1},'_Alpha_InducedPow_ITI.mat']);

Pre_SMD_Beta    = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{1},'_Beta_InducedPow_ITI.mat']);
Pos_SMD_Beta    = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{1},'_Beta_InducedPow_ITI.mat']);

Pre_SMD_GammaLo = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{1},'_GammaLo_InducedPow_ITI.mat']);
Pos_SMD_GammaLo = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{1},'_GammaLo_InducedPow_ITI.mat']);

Pre_SMD_GammaHi = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{1},'_GammaHi_InducedPow_ITI.mat']);
Pos_SMD_GammaHi = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{1},'_GammaHi_InducedPow_ITI.mat']);

subplot(4,3,2); hold on;
ForLineFigure(Pre_SMD_Alpha(:,8:end),T,'c'); 
ForLineFigure(Pos_SMD_Alpha(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;ylim([-1 16]); 
title('MD-8-14 (Hz)');hold off;
clear Pre_SMD_Alpha Pos_SMD_Alpha;


subplot(4,3,5); hold on;
ForLineFigure(Pre_SMD_Beta(:,8:end),T,'c'); 
ForLineFigure(Pos_SMD_Beta(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;ylim([-1 10]); 
title('MD-16-32 (Hz)');hold off;
clear Pre_SMD_Beta Pos_SMD_Beta;


subplot(4,3,8); hold on
ForLineFigure(Pre_SMD_GammaLo(:,8:end),T,'c'); 
ForLineFigure(Pos_SMD_GammaLo(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; ylim([-2.5 2.5]);
title('MD-33-64 (Hz)');hold off;
clear Pre_SMD_GammaLo Pos_SMD_GammaLo;

subplot(4,3,11); hold on
ForLineFigure(Pre_SMD_GammaHi(:,8:end),T,'c'); 
ForLineFigure(Pos_SMD_GammaHi(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)');axis tight; ylim([-1 2.2]);
title('MD-65-100 (Hz)');hold off;
clear Pre_SMD_GammaHi Pos_SMD_GammaHi;

%% Posterior-Ventral
Pre_SPV_Alpha   = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{3},'_Alpha_InducedPow_ITI.mat']);
Pos_SPV_Alpha   = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{3},'_Alpha_InducedPow_ITI.mat']);

Pre_SPV_Beta    = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{3},'_Beta_InducedPow_ITI.mat']);
Pos_SPV_Beta    = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{3},'_Beta_InducedPow_ITI.mat']);

Pre_SPV_GammaLo = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{3},'_GammaLo_InducedPow_ITI.mat']);
Pos_SPV_GammaLo = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{3},'_GammaLo_InducedPow_ITI.mat']);

Pre_SPV_GammaHi = importdata(['ForLinePlots/',Task{2},'_',Training{1},'_',Area{3},'_GammaHi_InducedPow_ITI.mat']);
Pos_SPV_GammaHi = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{3},'_GammaHi_InducedPow_ITI.mat']);

subplot(4,3,3); hold on;
ForLineFigure(Pre_SPV_Alpha(:,8:end),T,'c'); 
ForLineFigure(Pos_SPV_Alpha(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;ylim([-1 16]); 
 title('PV-8-14 (Hz)');hold off;
clear Pre_SPV_Alpha Pos_SPV_Alpha;

subplot(4,3,6); hold on;
ForLineFigure(Pre_SPV_Beta(:,8:end),T,'c'); 
ForLineFigure(Pos_SPV_Beta(:,8:end),T,'m');

xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;ylim([-1 10]); 
title('PV-16-32 (Hz)');hold off;
clear Pre_SPV_Beta Pos_SPV_Beta;


subplot(4,3,9); hold on
ForLineFigure(Pre_SPV_GammaLo(:,8:end),T,'c'); 
ForLineFigure(Pos_SPV_GammaLo(:,8:end),T,'m');

xline(-1,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; ylim([-2.5 2.5]);
title('PV-33-64 (Hz)');hold off;
clear Pre_SPV_GammaLo Pos_SPV_GammaLo;

subplot(4,3,12); hold on
ForLineFigure(Pre_SPV_GammaHi(:,8:end),T,'c'); 
ForLineFigure(Pos_SPV_GammaHi(:,8:end),T,'m');

xline(-1,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)');axis tight; ylim([-1 2.2]);
title('PV-65-100 (Hz)');hold off;
clear Pre_SPD_GammaHi Pos_SPD_GammaHi;

set(gcf,'renderer','Painters');set(gcf,'name',['FeatureTask_Induced']);
% print(fig,['FigS2.eps'],'-depsc');
% saveas(fig,['FigS2.fig'],'fig' );
% % close all;  
