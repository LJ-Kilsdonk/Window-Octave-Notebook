function [N_end,N_overtime]=multi_migration_run(phenotype_distribution,imm_freq,time_lengthrun,invasion_stepsize, MultiM_opt)
%MULTI_MIGRATION_RUN   creates a run with many random migration event,
%   the phenotype of each immigrant and the arrival time that is randomly
%   determined. The first arrival is at t=0. 
%   
%   [N_end,N_overtime]=MULTI_MIGRATION_RUN(phenotype_distribution,imm_freq,time_lengthrun,invasion_stepsize, MultiM_opt)
%   where N_end is a vector with the final population density of each
%   lineage, N_overtime is a matrix of population density per lineage
%   (columns) over time (rows), phenotype_distribution is a vector of the
%   probability distribution of phenotypes among the propagules,
%   imm_freq is the immigration frequency (i.e. the probability of an
%   immigration event per time step), time_lengthrun is a scalar with the
%   duration of the run, and MultiM_opt is a structure containing a list of
%   indeces of lineages under MultiM_opt.lin_j, the propagule density under
%   MultiM_opt.propsize and optionally information that helps to identify
%   the iteration in the loop from which the function is called under
%   MultiM_opt.indIter. 
%   
%   NOTE the probability that more than 1 propagule should enter during one
%   invasion_stepsize (a possibility ignored by the model) can be found
%   with: poisscdf(1,invasion_stepsize*imm_freq,'upper'). The expected
%   number of invasions in which more than 1 propagule should have been
%   entered during the run is thus:
%   poisscdf(1,inv_stepsize*disprate,'upper')*timeLrun  
% 
%   See also: updateQandd, runMgrow, loopmulti, preequilloopmulti

% 
global my_N my_t t a N n1 mut v_Ph a_orig MY_SETTINGS
% Simulations are faster if N starts with less columns. Growing
% N after each immigration event takes a bit of time, but numerically
% solving the model with a large matrix takes a lot of time.
n1=1; % number of lineages at start
t=0;
lin_j = MultiM_opt.lin_j; % indeces for invasion
propdensity = MultiM_opt.propsize;
n1max = length(lin_j); % maximum number of lineages
% total number of time indeces in output (including the first index t=0):
ntTotspan = time_lengthrun/invasion_stepsize +1; 
sim_time = linspace(0, time_lengthrun, ntTotspan); %time points in run
% use binomial distribution to determine arrival time of each migrant.
arrivalTime = find(binornd(1,imm_freq.*invasion_stepsize,...
    round(time_lengthrun./invasion_stepsize),1)).*invasion_stepsize;
% Set the first arrival at t = 0.
arrivalTime = [0;arrivalTime]; 
% randomly draw the phenotype for each migrant
[pheno_arriver,~] = find(mnrnd(1,phenotype_distribution,length(arrivalTime))');
% note: mnrnd is like a bernouili experiment, but with more than 2 possible
% outcomes. So instead of one p value it needs a vector of p values that
% add up to one. These p values are the phenotype_distribution values.

% farrM is a matrix of 2 collumns, first is the arrival time 2nd is the
% phenotype of arriver 
farrM=[arrivalTime, pheno_arriver]; % Founder arrival matrix


% INTRODUCE FIRST MIGRANT AND CORRECT DIMENSION N
%
% Without mutations only one phenotype is possible per lineage. This fact
% is utilized to speed up simulations without mutations, i.e. a vectorized
% version of N, where N is a column vector containing each lineage, and
% v_Ph is a column vector containg the phenotype of each lineage. The
% function updateQandd uses v_Ph to turn d into a column vector with as
% elements the resource-independent death rate of each lineage. As the
% calculation of mismatch depends on the number of phenotypes, the number
% of phenotypes has to be provided to updateQandd. This no longer
% corresponds to 'a', which is already the number of rows in N (i.e. 1).
% Therefore, another variable is used: a_orig to store this value.
if mut ==0
    % store number of phenotypes
    a_orig = a;
    a = 1;    
    % Set density of first immigrant to propsize.
    % Dimension of N will be 1 x 1
    N = propdensity; 
    
else
    firstmigrant = pheno_arriver(1); % phenotype of the first migrant
    % Set density of first immigrant to propsize.
    % Dimension of N will be a x n1 = a x 1
    enterfounder(firstmigrant,propdensity); 
end
MY_SETTINGS.Nelements = numel(N);

% Remove any immigrantion events that occur after the length of the
% simulation (this includes an invation exactly at t= time_lengthrun.
farrM( farrM(:,1) >= time_lengthrun, :) = [];

% Some dimensions for later:
n_inv = size(farrM,1); % number of invasions events (i.e. immigration events excluding that of the first migrant)

% if the maximum number of invasion is too much, display a warning and
% remove immigration events.
if n_inv > n1max
    warning('%d invasions but only %d possible by model, now reducing invasion_events ...',n_inv,n1max)
    farrM(n1max+1:end,:)=[];
    n_inv = size(farrM,1);
    if isfield(MultiM_opt,'indIter')
        disp(MultiM_opt.indIter);
    end
end

% Set phenotypes of each lineage in vector-N version
if mut ==0 % if vector N is used.
    % save all the phenotypes (the second column of farrM) to v_Ph,
    v_Ph = farrM(:,2); 
    % Later the second column of farrM is used to index the row in N, where
    % an immigrant propagule has to be added. Because N is now a column
    % vector, this row-index should always be 1.
    farrM(:,2) = 1;
end


% Allocate memory for output and save initial values
my_Ncomb = zeros(ntTotspan, n1max*(a));    
my_Ncomb(1,1:numel(N)) = N(:); 

%----------MAIN LOOP------------
i_tst  = 1; % time index of previous immigration event
if n_inv > 1
    for i_inv = 2:n_inv        
        % to avoid rounding errors, we here look for closest time point in
        % output time, to the event time, instead of time point that is
        % exactly identical.
        [~,i_tend] = min(abs( sim_time - farrM(i_inv,1) )); % time index of the next invasion
        
        % We want the time points tspan=sim_time(i_t0:i_tend), but if
        % length(tspan)==2, instead of just 2 time points, it gives all
        % internally calculated time points between these points. Therefore,
        % if length(tspan)==2 we add an intermediate time point in tspan.
        % And remove it afterwards.
        if i_tend - i_tst == 1 %i.e. length(tspan)==2
            runMgrow(sim_time(i_tst) : invasion_stepsize/2 : sim_time(i_tend));
            my_t(2)=[];
            my_N(2,:)=[];
        else
            runMgrow(sim_time(i_tst:i_tend));
        end
        
        
        % The previous run ended at the same time point, where the current
        % run started, the only difference being that the current run
        % starts after the migrant was added. To prevent double time
        % points, we remove the first time point from the current time
        % section. Thus, in the output, N is measured before migrants are
        % added and the migrant only appears in the output the time step
        % following the invasion (with exception of the first immigrant;
        % for which the first time step was saved in my_Ncomb before the
        % loop).  
        my_N(1,:) = []; 
        
        % add current run to combined run
        my_Ncomb(i_tst +1 : i_tend,1:size(my_N,2)) = my_N; 
        
        % last index of current run becomes the 1st of the next run.
        i_tst = i_tend; 
        
        %Do immigration event:
        N(farrM(i_inv,2),lin_j(i_inv)) = N(farrM(i_inv,2),lin_j(i_inv)) + propdensity;
    end
       
    % Provided the last immigration event didn't happen at the end of time_lengthrun,
    % run the till time_lengtrun (=sim_time(end))
    if i_tend < ntTotspan
        if ntTotspan - i_tst == 1 % i.e. length(tspan)==2
            runM(sim_time(i_tst) : invasion_stepsize/2 : sim_time(end),true);
            my_N(2,:)=[];
        else
            runM(sim_time(i_tst:end),true);
        end
        my_N(1,:) = [];
        my_Ncomb(i_tst +1 : end,1:size(my_N,2)) = my_N;
    end
    % use my_N to refer to N during the entire run
    my_N = my_Ncomb;
else % if no invasion
    runM(sim_time, true);
end
% use my_t to refer t during the entire run
my_t = sim_time'; 

% restore to the normal matrix version of N
if mut ==0 % if vector-N was used
    a = a_orig;
    % give my_Ncomb correct size
    my_Ncomb = zeros(ntTotspan, n1max*(a));
    % index of first element minus 1 of each column in matrix N (indexing
    % the matrix with one index, as is done in my_N)
    zerothelement = (0:n_inv-1)*a_orig;
    % index to turn  my_N(t,j) back to a full size (a x n1)-matrix my_N(t,ij)
    i_forN = zerothelement  + v_Ph'; 
    % turn my_N for vector-N to my_N as it would have been with matrix-N
    my_Ncomb(:,i_forN) = my_N(:,1:n_inv); 
    my_N = my_Ncomb;
end


% Process output
% Store last value of my_N to N. Because runMgrow.m does not add an extra
% column to my_N, any extra column it added to N will here be removed.
store

% Create matrix of final pop. density N(i,j), with i: phenotypes, j:
% lineage (group)
N_final = N; 
% Create vector of final pop density per lineage: N(j)
N_end=sum(N_final, 1); 
if nargout > 1
    % Create array of pop. densities over time N_array(t,i,j), with t:
    % time, i: phenotype, j: lineage (group)
    N_array = reshape(my_N, size(my_N,1),(a),size(my_N,2)/(a)); 
    % Create matrix of pop. density per lineage over time: N_overtime(t,j)
    N_overtime=squeeze(sum(N_array,2)); 
end

