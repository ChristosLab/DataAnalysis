%%%%%%%%%%%plot prefer location change 
load PPC_gauss_cue.mat
PPC_cue_pref=store_gauss_cue;
load PPC_gauss_delay.mat
PPC_delay_pref=store_gauss_delay;
load PFC_gauss_cue.mat
PFC_cue_pref=store_gauss_cue;
load PFC_gauss_delay.mat
f=figure;
PFC_delay_pref=store_gauss_delay;
scatter(PFC_cue_pref,PFC_delay_pref,60,[0 0.447 0.741],'LineWidth',1.5); 
hold on;
scatter(PPC_cue_pref,PPC_delay_pref,60,[0.8500 0.3250 0.0980],'LineWidth',1.5); 
xlim([-2,2]);
ylim([-2,2]);
xticklabels({'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'});
yticklabels({'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'});
hold on;
plot([-2,2],[-2,2],'k','LineWidth',1);
hAx=gca; 
hAx.LineWidth=1.5;hAx.FontSize = 12;
PFC_distanc=point_to_line_vec(1,-1,0,PFC_cue_pref,PFC_delay_pref);
PPC_distanc=point_to_line_vec(1,-1,0,PPC_cue_pref,PPC_delay_pref);
f.Position = [200 200 450 400];
lgd=legend('PFC','PPC');
lgd.FontSize = 14;
figure;
histogram(PFC_distanc,[-3:0.2:3],'Normalization','probability');
hold on;
histogram(PPC_distanc,[-3:0.2:3],'Normalization','probability');
 set(gca,'visible','off')