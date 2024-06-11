%% This script for Figure 3A_B
clc;  clear all; close all;
default_dir = pwd;
% 
Task     = {'Spatial','Feature'};
Training = {'PRETRAINING','POSTTRAINING'};
Area     = {'Mid-Dorsal','Posterior-Dorsal','Posterior-Ventral'};



Time = importdata(['ForLinePlots/TimeWhole.mat']);
T    = linspace(-2,5,length(Time));
T    = T(10:end);
%% Spatial
SPD_Alpha   = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{2},'_Alpha_InducedPow_ITI.mat']);
SPD_Beta    = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{2},'_Beta_InducedPow_ITI.mat']);
SPD_GammaLo = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{2},'_GammaLo_InducedPow_ITI.mat']);
SPD_GammaHi = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{2},'_GammaHi_InducedPow_ITI.mat']);
                        
SMD_Alpha   = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{1},'_Alpha_InducedPow_ITI.mat']);
SMD_Beta    = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{1},'_Beta_InducedPow_ITI.mat']);                           
SMD_GammaLo = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{1},'_GammaLo_InducedPow_ITI.mat']);                            
SMD_GammaHi = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{1},'_GammaHi_InducedPow_ITI.mat']);   
                      
SPV_Alpha   = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{3},'_Alpha_InducedPow_ITI.mat']);
SPV_Beta    = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{3},'_Beta_InducedPow_ITI.mat']);
SPV_GammaLo = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{3},'_GammaLo_InducedPow_ITI.mat']);  
SPV_GammaHi = importdata(['ForLinePlots/',Task{1},'_',Training{2},'_',Area{3},'_GammaHi_InducedPow_ITI.mat']);                            



fig = figure;
set(fig, 'Position', [30 40 850 650]);



subplot(4,2,1); hold on;
ForLineFigure(SPD_Alpha(:,10:end),T,'b'); 
ForLineFigure(SMD_Alpha(:,10:end),T,'r'); 
ForLineFigure(SPV_Alpha(:,10:end),T,'g'); 

% xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;%ylim([-15 15]); 
 title('PD-8-14 (Hz)');hold off;
clear SPD_Alpha SMD_Alpha SPV_Alpha;
% % 
subplot(4,2,3); hold on;
ForLineFigure(SPD_Beta(:,10:end),T,'b'); 
ForLineFigure(SMD_Beta(:,10:end),T,'r');
ForLineFigure(SPV_Beta(:,10:end),T,'g');

% xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;%ylim([-15 15]); 
title('PD-16-32 (Hz)');hold off;
clear SPD_Beta SMD_Beta SPV_Beta;
% % % 
subplot(4,2,5); hold on
ForLineFigure(SPD_GammaLo(:,10:end),T,'b'); 
ForLineFigure(SMD_GammaLo(:,10:end),T,'r'); 
ForLineFigure(SPV_GammaLo(:,10:end),T,'g'); 

% xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; %ylim([-6 1]);
title('PD-33-64 (Hz)');hold off;
clear SPD_GammaLo SMD_GammaLo SPV_GammaLo;

subplot(4,2,7); hold on
ForLineFigure(SPD_GammaHi(:,10:end),T,'b'); 
ForLineFigure(SMD_GammaHi(:,10:end),T,'r'); 
ForLineFigure(SPV_GammaHi(:,10:end),T,'g'); 

% xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)');axis tight; %ylim([-6 1]);
title('PD-65-100 (Hz)');hold off;
clear SPD_GammaHi SMD_GammaHi SPV_GammaHi;

%% Feature

FPD_Alpha   = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{2},'_Alpha_InducedPow_ITI.mat']);
FPD_Beta    = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{2},'_Beta_InducedPow_ITI.mat']);    
FPD_GammaLo = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{2},'_GammaLo_InducedPow_ITI.mat']);
FPD_GammaHi = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{2},'_GammaHi_InducedPow_ITI.mat']);
                      
FMD_Alpha   = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{1},'_Alpha_InducedPow_ITI.mat']);
FMD_Beta    = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{1},'_Beta_InducedPow_ITI.mat']);                          
FMD_GammaLo = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{1},'_GammaLo_InducedPow_ITI.mat']);                
FMD_GammaHi = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{1},'_GammaHi_InducedPow_ITI.mat']);

FPV_Alpha   = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{3},'_Alpha_InducedPow_ITI.mat']);
FPV_Beta    = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{3},'_Beta_InducedPow_ITI.mat']);
FPV_GammaLo = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{3},'_GammaLo_InducedPow_ITI.mat']);
FPV_GammaHi = importdata(['ForLinePlots/',Task{2},'_',Training{2},'_',Area{3},'_GammaHi_InducedPow_ITI.mat']);                           
subplot(4,2,2); hold on;
ForLineFigure(FPD_Alpha(:,10:end),T,'b'); 
ForLineFigure(FMD_Alpha(:,10:end),T,'r'); 
ForLineFigure(FPV_Alpha(:,10:end),T,'g'); 

% xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;%ylim([-15 15]); 
 title('PD-8-14 (Hz)');hold off;
clear FPD_Alpha FMD_Alpha FPV_Alpha;
% % 
subplot(4,2,4); hold on;
ForLineFigure(FPD_Beta(:,10:end),T,'b'); 
ForLineFigure(FMD_Beta(:,10:end),T,'r');
ForLineFigure(FPV_Beta(:,10:end),T,'g');

% xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight;%ylim([-15 15]); 
title('PD-16-32 (Hz)');hold off;
clear FPD_Beta FMD_Beta FPV_Beta;
% % % 
subplot(4,2,6); hold on
ForLineFigure(FPD_GammaLo(:,10:end),T,'b'); 
ForLineFigure(FMD_GammaLo(:,10:end),T,'r'); 
ForLineFigure(FPV_GammaLo(:,10:end),T,'g'); 

% xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)'); axis tight; %ylim([-6 1]);
title('PD-33-64 (Hz)');hold off;
clear SPD_GammaLo SMD_GammaLo SPV_GammaLo;

subplot(4,2,8); hold on
ForLineFigure(FPD_GammaHi(:,10:end),T,'b'); 
ForLineFigure(FMD_GammaHi(:,10:end),T,'r'); 
ForLineFigure(FPV_GammaHi(:,10:end),T,'g'); 

% xline(-1, ':','LineWidth',1,'Color',[0 0 0]);
xline(0,  ':','LineWidth',1,'Color',[0 0 0]);
xline(0.5,':','LineWidth',1,'Color',[0 0 0]);
xline(2.0,':','LineWidth',1,'Color',[0 0 0]);
xline(2.5,':','LineWidth',1,'Color',[0 0 0]);
xline(4.0,':','LineWidth',1,'Color',[0 0 0]);
xlabel('Time'); ylabel('Power(dB)');axis tight; %ylim([-6 1]);
title('PD-65-100 (Hz)');hold off;
clear FPD_GammaHi FMD_GammaHi FPV_GammaHi;
 
set(gcf,'renderer','Painters');set(gcf,'name',['Spatial_FeatureTask_Areas_Induced']);
% print(fig,['FigS4.eps'],'-depsc');
% saveas(fig,['FigS4.fig'],'fig' );
% close all;  
