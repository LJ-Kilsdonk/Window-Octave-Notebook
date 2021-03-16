function [n1max,distrn1max]=n1maxtest(imm_freq,time_lengthrun,invasion_stepsize, m_repl)
%N1MAXTEST  Tests the maximum number of lineages in a simulation.
%   In the simulations m_repl replicates are used. To see the expected
%   distribution of the maximum number of lineages, the functions repeats
%   the calculation of arrival times of immigrants in each of the m_repl
%   replicates 1000 times. It plots this distribution. 
%   n1max contains the overall maximum lineage richness. distrn1max
%   contains the distribution of n1max values among the 1000 replicates of
%   m_repl simulations.

reps  = 1e3; % number of times to repeat 

distrn1max=zeros(reps,1); % distribution of maximum number of lineages
for i=1:reps
    n1m=0;
    for j=1:m_repl
        n1 = length(find(binornd(1,imm_freq.*invasion_stepsize,round(time_lengthrun./invasion_stepsize),1)).*invasion_stepsize);
        n1m = max(n1m, n1+1);
    end
    distrn1max(i)=n1m;
    
end
n1max=max(distrn1max);

% plot distribution of n1max.
figure
histogram(distrn1max,'Normalization','probability')
xlabel('n1max')
ylabel('probability')
