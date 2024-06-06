function [s, m, sem, loctt, rateyy, d, R] = gaus_fit_8_loc(input_matrix)
% Fits 1-D Gaussian curve to 8 spatial locations
% Input format: n by 9 with peak at 5
% based on Christos Constantinidis, 11-Mar-2010
% Junda Zhu, 20230620

loc=1:9;
if min(size(input_matrix)) >= 1
    [s, m, sem, ~] = allmeans(1:9,input_matrix);
    % rate(1:9)=inp_rate(1:9);
    rate(1:9) = m(1:9);

    % 09-Aug-2011, xzhou
    myfunc=inline('beta(1)+beta(2)*exp(-0.5 *((loc - beta(3))/beta(4)) .^2)','beta','loc');

    beta0 = [min(rate) max(rate)-min(rate) 5 1];

    beta=nlinfit(loc,rate,myfunc,beta0);
    a=beta(1);b=beta(2);c=beta(3);
    d=beta(4);

    loctt=min(loc):0.1:max(loc);
    rateyy=a+b*exp(-0.5 *((loctt - c)/d) .^2);

    % calculate R (corr)
    loc = 1:9;
    rateyy2=a+b*exp(-0.5 *((loc - c)/d) .^2);
    R = corr(rate(1:9)', rateyy2');
else
    [s, m, sem, loctt, rateyy, d, R] = deal(nan(1,9));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s, m, sem, nt] = allmeans(st, r)
% Given the stimulus (st) and the responses (r) at each trial -- trials
% correspond to rows in r -- this function returns a vector of stimuli
% without repeats and means, SEMs and numbers of trials.
st = st(:);
% try
%     ii = find(isnan(st)==1 | isnan(r)==1);
%     st(ii) = [];
%     r(ii)  = [];
% end
[s, ~] = unique(st);
Ns = length(s);
for j=1:Ns
    ii = find(~isnan(r(j,:)));
    m(j)  = nanmean(r(j,:));
    sem(j) = std(r(j,:),'omitnan')./sqrt(length(ii));     % this can be zero
    nt1(j) = length(ii);
end
nt = min(nt1);
m = m(:);
sem = sem(:);