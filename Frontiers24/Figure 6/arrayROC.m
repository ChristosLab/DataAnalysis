function [areaROC, prob, criter] = arrayROC(classtarget,classdistr)

% Constructs ROC curve and computes area under it
% Compares two sets of firing rates (target and distractor) and
% determines how discriminable they are relative to each other
%
% Author: Christos Constantinidis, Ph.D.
% 11-MAR-2009

hits=length(classtarget);
fals=length(classdistr);

allrates=[classtarget classdistr];
criter=unique(allrates);
tot=length(criter);

for i=1:tot
    hitprob=sum(classtarget>=criter(i));
    falsprob=sum(classdistr>=criter(i));
    prob(i,1)=hitprob./hits;
    prob(i,2)=falsprob./fals;
end
% ROC curve includes 0,0 point
prob(tot+1,1)=0;
prob(tot+1,2)=0;

% plot(prob(:,2), prob(:,1), '-o');

areaROC=0;
for i=1:tot
    areaROC=areaROC + (prob(i,1)+prob(i+1,1))*(prob(i,2)-prob(i+1,2))/2;
end


