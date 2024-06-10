%% This script for relation gamma power and spike rate
clc; clear;close all;
default_dir=pwd;


Task     = {'Spatial','Feature'};
Training = {'PRETRAINING','POSTTRAINING'};
INF      = {'SELECTIVE','NON_SELECTIVE','PARTIAL_SELECTIVE'};

bin_edges = [-1.5 -1 -0.7 -0.4 -0.1 0.2 0.5 0.8 1.1 1.4 1.7 2 2.3]; %
for jj=1:length(bin_edges)-1
    xx(jj)=(bin_edges(jj)+bin_edges(jj+1))/2;
end


sz = 20;

TEMP_nonsel = readtable('Data/LogRFpower/Dataforplot_Spatial_pretraining_firingrate_Epoch_Nsl.csv');
TEMP_par    = readtable('Data/LogRFpower/Dataforplot_Spatial_pretraining_firingrate_Epoch_Par.csv');
TEMP_sel    = readtable('Data/LogRFpower/Dataforplot_Spatial_pretraining_firingrate_Epoch_Sel.csv');


subplot(2,2,1);
scatter(TEMP_nonsel.Firingrate,TEMP_nonsel.GammaPow,'MarkerEdgeColor',[0.6 0.6 1]); hold on;   %%%% 'MarkerFaceColor',[0.6 0.6 1]
scatter(TEMP_par.Firingrate,TEMP_par.GammaPow,'MarkerEdgeColor',[0.7 0.7 0.7]);
scatter(TEMP_sel.Firingrate,TEMP_sel.GammaPow,'MarkerEdgeColor',[1 0.9 1]);



[Nnsl,Edgensl,binnsl]     = histcounts(TEMP_nonsel.Firingrate,bin_edges);
[Gmeannslval, Gstdnslval] = grpstats(TEMP_nonsel.GammaPow,binnsl,{@mean,@std});
xx_nsl                    = xx(find(Nnsl>0));
Ind                       = find(Gstdnslval>0);
scatter(xx_nsl(Ind),Gmeannslval(Ind),sz,'filled','b');
errorbar(xx_nsl(Ind), Gmeannslval(Ind),Gstdnslval(Ind), 'color','b','LineStyle','none'); 
Myrefline = refline([TEMP_nonsel.Firingrate_1(1),TEMP_nonsel.X_Intercept_(1)]);
Myrefline.Color     = 'b';
Myrefline.LineWidth = 1;
Myrefline.XData     = [-2 3];
clear TEMP_nonsel Nnsl Edgensl binnsl Gmeannslval Gstdnslval xx_nsl Ind Myrefline;
%  
 
 
[Npar,Edgepar,binpar]     = histcounts(TEMP_par.Firingrate,bin_edges);
[Gmeanparval, Gstdparval] = grpstats(TEMP_par.GammaPow,binpar,{@mean,@std});
xx_par                    = xx(find(Npar>0));
Ind                       = find(Gstdparval>0);
scatter(xx_par(Ind),Gmeanparval(Ind),sz,'filled','k'); 
errorbar(xx_par(Ind), Gmeanparval(Ind),Gstdparval(Ind), 'color','k','LineStyle','none');
Myrefline           = refline([TEMP_par.Firingrate_1(1),TEMP_par.X_Intercept_(1)]);
Myrefline.Color     = 'k';
Myrefline.LineWidth = 1; 
Myrefline.XData     = [-2 3];
clear TEMP_par Npar Edgepar binpar Gmeanparval Gstdparval xx_par Myrefline Ind;


[Nsel,Edgesel,binsel]     = histcounts(TEMP_sel.Firingrate,bin_edges);
[Gmeanselval, Gstdselval] = grpstats(TEMP_sel.GammaPow,binsel,{@mean,@std});
xx_sel                    = xx(find(Nsel>0));
Ind                       = find(Gstdselval>0);

scatter(xx_sel(Ind),Gmeanselval(Ind),sz,'filled','m'); hold on
errorbar(xx_sel(Ind), Gmeanselval(Ind),Gstdselval(Ind),'color','m', 'LineStyle','none'); 
Myrefline = refline([TEMP_sel.Firingrate_1(1),TEMP_sel.X_Intercept_(1)]);
Myrefline.Color = 'm';
Myrefline.LineWidth = 1; 
Myrefline.XData     = [-2 3];
clear TEMP_sel Nsel Edgesel binsel Gmeanselval Gstdselval xx_sel Ind Myrefline; 



% legend('NonInfo','','','Partially','','','Info');
title([Task{1},'-',Training{1}]);
xlabel('log(Firing rate [Sp/S])'); ylabel('Raw gamma power(dB)');
xlim([-2 3]); ylim([10 55]);
clear TEMP_sel TEMP_par TEMP_nonsel; 
hold off;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


TEMP_nonsel  = readtable('Data/LogRFpower/Dataforplot_Spatial_pretraining_firingrate_Epochwithbaseline_Nsl.csv');
TEMP_par     = readtable('Data/LogRFpower/Dataforplot_Spatial_pretraining_firingrate_Epochwithbaseline_Par.csv');
TEMP_sel     = readtable('Data/LogRFpower/Dataforplot_Spatial_pretraining_firingrate_Epochwithbaseline_Sel.csv');


subplot(2,2,3);

scatter(TEMP_nonsel.Firingrate,TEMP_nonsel.GammaPow,'MarkerEdgeColor',[0.6 0.6 1]); hold on;
scatter(TEMP_par.Firingrate,TEMP_par.GammaPow,'MarkerEdgeColor',[0.7 0.7 0.7]);
scatter(TEMP_sel.Firingrate,TEMP_sel.GammaPow,'MarkerEdgeColor',[1 0.9 1]);




[Nnsl,Edgensl,binnsl]     = histcounts(TEMP_nonsel.Firingrate,bin_edges);
[Gmeannslval, Gstdnslval] = grpstats(TEMP_nonsel.GammaPow,binnsl,{@mean,@std});
xx_nsl                    = xx(find(Nnsl>0));
Ind                       = find(Gstdnslval>0);
scatter(xx_nsl(Ind),Gmeannslval(Ind),sz,'filled','b'); 
errorbar(xx_nsl(Ind), Gmeannslval(Ind),Gstdnslval(Ind), 'color','b','LineStyle','none'); 
Myrefline = refline([TEMP_nonsel.Firingrate_1(1),TEMP_nonsel.X_Intercept_(1)]);
Myrefline.Color     = 'b';
Myrefline.LineWidth = 1;
Myrefline.XData     = [-0.5 2.5];
clear TEMP_nonsel Nnsl Edgensl binnsl Gmeannslval Gstdnslval xx_nsl Ind Myrefline;



[Npar,Edgepar,binpar]     = histcounts(TEMP_par.Firingrate,bin_edges);
[Gmeanparval, Gstdparval] = grpstats(TEMP_par.GammaPow,binpar,{@mean,@std});
xx_par                    = xx(find(Npar>0));
Ind                       = find(Gstdparval>0);

scatter(xx_par(Ind),Gmeanparval(Ind),sz,'filled','k');
errorbar(xx_par(Ind), Gmeanparval(Ind),Gstdparval(Ind), 'color','k','LineStyle','none');
Myrefline           = refline([TEMP_par.Firingrate_1(1),TEMP_par.X_Intercept_(1)]);
Myrefline.Color     = 'k';
Myrefline.LineWidth = 1; 
Myrefline.XData     = [-0.5 2.5];
clear TEMP_par Npar Edgepar binpar Gmeanparval Gstdparval xx_par Myrefline Ind;


[Nsel,Edgesel,binsel]     = histcounts(TEMP_sel.Firingrate,bin_edges);
[Gmeanselval, Gstdselval] = grpstats(TEMP_sel.GammaPow,binsel,{@mean,@std});
xx_sel                    = xx(find(Nsel>0));
Ind                       = find(Gstdselval>0);

scatter(xx_sel(Ind),Gmeanselval(Ind),sz,'filled','m');
errorbar(xx_sel(Ind), Gmeanselval(Ind),Gstdselval(Ind),'color','m', 'LineStyle','none'); 
Myrefline = refline([TEMP_sel.Firingrate_1(1),TEMP_sel.X_Intercept_(1)]);
Myrefline.Color = 'm';
Myrefline.LineWidth = 1; 
Myrefline.XData     = [-0.5 2.5];
clear TEMP_sel Nsel Edgesel binsel Gmeanselval Gstdselval xx_sel Ind Myrefline; 




% legend('NonInfo','','','Partially','','','Info');
title([Task{1},'-',Training{1}]);
xlabel('log(Firing rate [Sp/S])'); ylabel('Gamma power-Baseline (dB)');
xlim([-0.5 2.5]);ylim([-3 10]);

clear TEMP_sel TEMP_par TEMP_nonsel; 
hold off;


%%





bin_edges = [-2 -1.5 -1 -0.7 -0.4 -0.1 0.2 0.5 0.8 1.1 1.4 1.7 2 2.5]; %
for jj=1:length(bin_edges)-1
    xx(jj)=(bin_edges(jj)+bin_edges(jj+1))/2;
end



TEMP_nonsel               = readtable('Data/LogRFpower/Dataforplot_Spatial_posttraining_firingrate_Epoch_Nsl.csv');
TEMP_par                  = readtable('Data/LogRFpower/Dataforplot_Spatial_posttraining_firingrate_Epoch_Par.csv');
TEMP_sel                  = readtable('Data/LogRFpower/Dataforplot_Spatial_posttraining_firingrate_Epoch_Sel.csv');


subplot(2,2,2);
scatter(TEMP_nonsel.Firingrate,TEMP_nonsel.GammaPow,'MarkerEdgeColor',[0.6 0.6 1]); hold on;
scatter(TEMP_par.Firingrate,TEMP_par.GammaPow,'MarkerEdgeColor',[0.7 0.7 0.7]);
scatter(TEMP_sel.Firingrate,TEMP_sel.GammaPow,'MarkerEdgeColor',[1 0.9 1]);



[Nnsl,Edgensl,binnsl]     = histcounts(TEMP_nonsel.Firingrate,bin_edges);
[Gmeannslval, Gstdnslval] = grpstats(TEMP_nonsel.GammaPow,binnsl,{@mean,@std});
xx_nsl                    = xx(find(Nnsl>0));
Ind                       = find(Gstdnslval>0);
scatter(xx_nsl(Ind),Gmeannslval(Ind),sz,'filled','b'); 
errorbar(xx_nsl(Ind), Gmeannslval(Ind),Gstdnslval(Ind), 'color','b','LineStyle','none'); 
Myrefline = refline([TEMP_nonsel.Firingrate_1(1),TEMP_nonsel.X_Intercept_(1)]);
Myrefline.Color     = 'b';
Myrefline.LineWidth = 1;
Myrefline.XData     = [-2 3];
clear  TEMP_nonsel Nnsl Edgensl binnsl Gmeannslval Gstdnslval xx_nsl Ind Myrefline;



[Npar,Edgepar,binpar]     = histcounts(TEMP_par.Firingrate,bin_edges);
[Gmeanparval, Gstdparval] = grpstats(TEMP_par.GammaPow,binpar,{@mean,@std});
xx_par                    = xx(find(Npar>0));
Ind                       = find(Gstdparval>0);

scatter(xx_par(Ind),Gmeanparval(Ind),sz,'filled','k');
errorbar(xx_par(Ind), Gmeanparval(Ind),Gstdparval(Ind), 'color','k','LineStyle','none');
Myrefline           = refline([TEMP_par.Firingrate_1(1),TEMP_par.X_Intercept_(1)]);
Myrefline.Color     = 'k';
Myrefline.LineWidth = 1; 
Myrefline.XData     = [-2 3];
clear TEMP_par Npar Edgepar binpar Gmeanparval Gstdparval xx_par Myrefline Ind;

[Nsel,Edgesel,binsel]     = histcounts(TEMP_sel.Firingrate,bin_edges);
[Gmeanselval, Gstdselval] = grpstats(TEMP_sel.GammaPow,binsel,{@mean,@std});
xx_sel                    = xx(find(Nsel>0));
Ind                       = find(Gstdselval>0);

scatter(xx_sel(Ind),Gmeanselval(Ind),sz,'filled','m');
errorbar(xx_sel(Ind), Gmeanselval(Ind),Gstdselval(Ind),'color','m', 'LineStyle','none'); 
Myrefline = refline([TEMP_sel.Firingrate_1(1),TEMP_sel.X_Intercept_(1)]);
Myrefline.Color = 'm';
Myrefline.LineWidth = 1; 
Myrefline.XData     = [-2 3];
clear TEMP_sel Nsel Edgesel binsel Gmeanselval Gstdselval xx_sel Ind Myrefline; 

title([Task{1},'-',Training{2}]);
xlabel('log(Firing rate [Sp/S])'); ylabel('Raw gamma power(dB)');
xlim([-2 3]); ylim([10 55]);
clear TEMP_sel TEMP_par TEMP_nonsel; 
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TEMP_nonsel               = readtable('Data/LogRFpower/Dataforplot_Spatial_posttraining_firingrate_Epochwithbaseline_Nsl.csv');
TEMP_par                  = readtable('Data/LogRFpower/Dataforplot_Spatial_posttraining_firingrate_Epochwithbaseline_Par.csv');
TEMP_sel                  = readtable('Data/LogRFpower/Dataforplot_Spatial_posttraining_firingrate_Epochwithbaseline_Sel.csv');

subplot(2,2,4);
scatter(TEMP_nonsel.Firingrate,TEMP_nonsel.GammaPow,'MarkerEdgeColor',[0.6 0.6 1]); hold on
scatter(TEMP_par.Firingrate,TEMP_par.GammaPow,'MarkerEdgeColor',[0.7 0.7 0.7]);
scatter(TEMP_sel.Firingrate,TEMP_sel.GammaPow,'MarkerEdgeColor',[1 0.9 1]);



[Nnsl,Edgensl,binnsl]     = histcounts(TEMP_nonsel.Firingrate,bin_edges);
[Gmeannslval, Gstdnslval] = grpstats(TEMP_nonsel.GammaPow,binnsl,{@mean,@std});
xx_nsl                    = xx(find(Nnsl>0));
Ind                       = find(Gstdnslval>0);
scatter(xx_nsl(Ind),Gmeannslval(Ind),sz,'filled','b'); 
errorbar(xx_nsl(Ind), Gmeannslval(Ind),Gstdnslval(Ind), 'color','b','LineStyle','none'); 
Myrefline = refline([TEMP_nonsel.Firingrate_1(1),TEMP_nonsel.X_Intercept_(1)]);
Myrefline.Color     = 'b';
Myrefline.LineWidth = 1;
Myrefline.XData     = [-0.5 2.5];
clear TEMP_nonsel Nnsl Edgensl binnsl Gmeannslval Gstdnslval xx_nsl Ind Myrefline;



[Npar,Edgepar,binpar]     = histcounts(TEMP_par.Firingrate,bin_edges);
[Gmeanparval, Gstdparval] = grpstats(TEMP_par.GammaPow,binpar,{@mean,@std});
xx_par                    = xx(find(Npar>0));
Ind                       = find(Gstdparval>0);
scatter(xx_par(Ind),Gmeanparval(Ind),sz,'filled','k');
errorbar(xx_par(Ind), Gmeanparval(Ind),Gstdparval(Ind), 'color','k','LineStyle','none');
Myrefline           = refline([TEMP_par.Firingrate_1(1),TEMP_par.X_Intercept_(1)]);
Myrefline.Color     = 'k';
Myrefline.LineWidth = 1; 
Myrefline.XData     = [-0.5 2.5];
clear TEMP_par Npar Edgepar binpar Gmeanparval Gstdparval xx_par Myrefline Ind;

[Nsel,Edgesel,binsel]     = histcounts(TEMP_sel.Firingrate,bin_edges);
[Gmeanselval, Gstdselval] = grpstats(TEMP_sel.GammaPow,binsel,{@mean,@std});
xx_sel                    = xx(find(Nsel>0));
Ind                       = find(Gstdselval>0);
scatter(xx_sel(Ind),Gmeanselval(Ind),sz,'filled','m'); 
errorbar(xx_sel(Ind), Gmeanselval(Ind),Gstdselval(Ind),'color','m', 'LineStyle','none'); 
Myrefline = refline([TEMP_sel.Firingrate_1(1),TEMP_sel.X_Intercept_(1)]);
Myrefline.Color = 'm';
Myrefline.LineWidth = 1; 
Myrefline.XData     = [-0.5 2.5];
clear TEMP_sel  Nsel Edgesel binsel Gmeanselval Gstdselval xx_sel Ind Myrefline; 


title([Task{1},'-',Training{2}]);
xlabel('log(Firing rate [Sp/S])'); ylabel('Gamma power-Baseline (dB)');
xlim([-0.5 2.5]);ylim([-3 10]);
clear TEMP_sel TEMP_par TEMP_nonsel; 
hold off;

