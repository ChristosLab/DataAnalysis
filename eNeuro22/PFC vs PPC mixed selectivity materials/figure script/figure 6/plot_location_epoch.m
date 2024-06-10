clear;
load pfc_slide_f.mat
pfc_f=store_interaction_f;
load ppc_slide_f.mat
ppc_f=store_interaction_f;
pfc_f(find(pfc_f>20))=NaN;
ppc_f(find(ppc_f>20))=NaN;
pfc_result=nanmean(pfc_f');
ppc_result=nanmean(ppc_f');
plot(pfc_result);
hold on;
plot(ppc_result);
temp_c=polyfit([1:length(pfc_result)],pfc_result,1);
pfc_slope=temp_c(1);
temp_c=polyfit([1:length(ppc_result)],ppc_result,1);
ppc_slope=temp_c(1);