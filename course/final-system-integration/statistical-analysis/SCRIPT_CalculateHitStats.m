%% SCRIPT_CalculateHitStats
% This script compares the probability of an exact number of hits to the
% probability of *at least* a specified number of hits.
%   Variables:
%       N     - total number of shots
%       pOne  - probability of *at least* one hit when taking N shots
%       pHit  - probability that a single shot hits the target
%       p(i)  - probability that exactly i-1 shots out of N will hit the
%               target
%       pX(i) - probability that *at least* i-1 shots out of N will hit the
%               target
%
%   D. Evangelista & M. Kutzer, USNA, 17Apr2020
clear all
close all
clc

%% Define student parameters
N = 8;         % N shots
pOne = 0.70;    % 95% probability of at least a single hit

%% Calculate the probability of each shot hitting
pHit = pOne2pHit(pOne,N);
fprintf('Probability of a single shot hitting:     %8.4f%%\n',pOne*100);
fprintf('Probability of an individual shot hitting:%8.4f%%\n',pHit*100);
fprintf('\n');
%% Calculate probabilities for a given number of hits
for n = 0:N
    p(n+1) = prob_nhits(n,N,pHit);
    fprintf('\tProbability of %2d of %2d shots hitting the target: %8.4f%%\n',n,N,p(n+1)*100);
end
fprintf('\t\tTOTAL: %8.4f%%\n',sum(p)*100);
fprintf('\n');
%% Calculate probabilities for *at least* a given number of hits
for n = 0:N
    pX(n+1) = sum(p( (n+1):end ));
    fprintf('\tProbability of at least %2d of %2d shots hitting the target: %8.4f%%\n',n,N,pX(n+1)*100);
end

%% plot results
fig = figure;
axs = axes('Parent',fig);
hold(axs,'on');
xlabel('Number of shots (n) that hit the target');
ylabel('Probability');
for n = 0:N
    pltCUM(n+1) = plot(axs,n*[1,1],100*pX(n+1)*[0,1],'r','LineWidth',12);
    pltSNG(n+1) = plot(axs,n*[1,1],100* p(n+1)*[0,1],'b','LineWidth',10);
end
xlim(axs,[-1,N+1]);
ylim(axs,[0,100]);

yticks = 0:10:100;
for i = 1:numel(yticks)
    yticklabels{i} = sprintf('%d%%',yticks(i));
end
set(axs,'xtick',0:N,'ytick',yticks,'yticklabel',yticklabels);

legend(axs,{sprintf('Minimum Number of Hits\n(i.e. at least n hits)'),...
    sprintf('Exact Number of Hits\n(i.e. exactly n hits)')});

saveas(fig,sprintf('ProbabilityAnalysis_%dshots_%dpercentOfsingleHit.png',N,round(pOne*100)),'png');