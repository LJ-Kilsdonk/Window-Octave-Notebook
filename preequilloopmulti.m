function [linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,MultiM_opt] = preequilloopmulti(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file,MultiM_opt)
%PREEQUILLOOPMULTI  Main loop for pre-equilibrium runs with multiple
%   migrations
%   
%   Loops over different settings and replicates, each iteration it calls
%   multi_migration_run twice (once with and once without mutations) to
%   perform a simulation with multiple immigration events for which it
%   stores the lineage richness and arrival rank over time
% 
%   [linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,MultiM_opt] = 
%   PREEQUILLOOPMULTI(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file,MultiM_opt) 
%   where linR_Array is the lineage richness with mutations, linR_ArrayM
%   the lineage richness without mutations, gmeanARank_Array the mean
%   arrival rank with mutations, gmeanARank_Array_nomut the mean arrival
%   rank without mutations. For each the dimensions are time x E x
%   Immigration frequency.
%   m_repl is the number of replicates, timeLrun the duration of a run,
%   list_imm_freq a list of the immigration frequencies to be used, listE a
%   list of the E values to be tested, isredo a boolean operator that
%   states if the simulations are redone, even if name_file already exists.
%   If false and the file exists, the output data is loaded from the file.
%   MultiM_opt is a structure used to set optional values different from
%   the default. Default maximum number of lineage is 180. To change this,
%   change MultiM_opt.n1max
%   
%   The minimal density for a lineage to be counted in lineage richness
%   is 105% the propagule density.
% 
%   See also: loopmulti, multi_migration_run
global E propsize mut a

% default behavior if the file to load from/save to doesn't exist, yet
% isredo is false, is to do the simulations anyway and save them under
% filename.
if isredo == 1 || ~(exist(name_file,'file') == 2)
    MultiM_opt.warnings =[]; % not used in current version; but initializes MultiM_opt if not provided.
    
    % create probability distribution of phenotype frequencies among migrants
    amin1 = a - 1; 
    S= (factorial(amin1)./(factorial(0:amin1).*factorial(amin1 - (0:amin1))))./(2^amin1); 
    MultiM_opt.S = S;
    
    % previous versions allowed for grouping several immigrants under the
    % same 'lineage group'. Here each migrant corresponds to one lineage.
    % Therefore the column index of each lineage is simply the arrival rank
    % of that lineage.
    MultiM_opt.n1max = getknownfield(MultiM_opt,'n1max',180); % max number of lineages
    MultiM_opt.lin_j=(1:MultiM_opt.n1max)'; % column index of each lineage, ordered by arrival rank
    
    % set default propagule density
    if isempty(propsize)
        propsize = 0.001;
    end
    MultiM_opt.propsize = getknownfield(MultiM_opt,'propsize',propsize);
    propsize = MultiM_opt.propsize;
    
    Thr= propsize.*1.05; % minimum density to be included in the lineage richness count
    
    MultiM_opt.inv_stepsize = getknownfield(MultiM_opt,'inv_stepsize',0.1);
    inv_stepsize = MultiM_opt.inv_stepsize;
    
    % total number of time indeces in output multinvHigh (including the
    % first index t=0):
    ntTotspan = timeLrun / inv_stepsize +1; 
    % To reduce data points, and increase processing speed, we take
    % steps of 1 generation between the observations, and discard
    % the data at time points in between. Note this does not
    % influence the results, only the resolution of the figures, since 1
    % data point per time step is plotted, instead of 10 (=1/inv_stepsize)
    MultiM_opt.Roundedtime = 0:timeLrun;
    iRoundedtime = 1 : 1 / inv_stepsize : ntTotspan; % index in output time corresponding to Roundedtime.
    n_Rtime = length(MultiM_opt.Roundedtime);
    
    % Allocate memory
    linR_Array               =   zeros(n_Rtime,length(listE),length(list_imm_freq));
    gmeanARank_Array         =   linR_Array;
    linR_Array_nomut         =   linR_Array;
    gmeanARank_Array_nomut   =   linR_Array;
    
    tic
    % each iteration tests with and without mutations, thus at the start of
    % each iteration the original mutation rate needs to be restored.
    mutold=mut; % the original mutation rate.
    % MAIN LOOP
    for i_E=1:length(listE)
        E = listE(i_E);
        for i_imm_freq = 1:length(list_imm_freq)
            imm_freq = list_imm_freq(i_imm_freq);
            for i_rep=1:m_repl
                % restore mutation rate
                mut = mutold;
                % store information unique to current loop (for debugging)
                MultiM_opt.indIter=[E imm_freq i_rep mut];
                % Do a run with mutations to find density per lineage
                % over time
                [~,N_All]=multi_migration_run(S,imm_freq,timeLrun,inv_stepsize, MultiM_opt);
                % Reduce resolution
                N_AllRT=N_All(iRoundedtime,:);
                % Add current lineage richness (sum(N_AllRT>Thr,2)) to that
                % of previous replicate runs with same richness. After the
                % loop it is divided by m_repl to get the mean over the
                % replicates.
                linR_Array(:,i_E,i_imm_freq)= linR_Array(:,i_E,i_imm_freq) + sum(N_AllRT>Thr,2);
                % Same as linR_Array, but for the mean arrival rank. The
                % mean arrival rank is the average arrival rank weighted
                % by the frequency of the lineage in the population.
                gmeanARank_Array(:,i_E,i_imm_freq)= gmeanARank_Array(:,i_E,i_imm_freq) + sum((N_AllRT.*repmat(1:size(N_AllRT,2),n_Rtime,1)),2)./sum(N_AllRT,2);
                
                % Same as above, but without mutations
                mut=0;
                MultiM_opt.indIter=[E imm_freq i_rep mut];
                [~,N_All]=multi_migration_run(S,imm_freq,timeLrun,inv_stepsize, MultiM_opt);
                N_AllRT=N_All(iRoundedtime,:);
                linR_Array_nomut(:,i_E,i_imm_freq) = linR_Array_nomut(:,i_E,i_imm_freq) + sum(N_AllRT>Thr,2);
                gmeanARank_Array_nomut(:,i_E,i_imm_freq) = gmeanARank_Array_nomut(:,i_E,i_imm_freq) + sum((N_AllRT.*repmat(1:size(N_AllRT,2),n_Rtime,1)),2)./sum(N_AllRT,2);
                
                % display repetition to monitor progress
                disp(i_rep);
            end
            % Divide by m_repl to get the mean, instead of the sum
            linR_Array(:,i_E,i_imm_freq) = linR_Array(:,i_E,i_imm_freq) / m_repl;
            gmeanARank_Array(:,i_E,i_imm_freq) = gmeanARank_Array(:,i_E,i_imm_freq) / m_repl;
            linR_Array_nomut(:,i_E,i_imm_freq) = linR_Array_nomut(:,i_E,i_imm_freq) / m_repl;
            gmeanARank_Array_nomut(:,i_E,i_imm_freq) = gmeanARank_Array_nomut(:,i_E,i_imm_freq) / m_repl;
            
            % display settings to monitor progress
            disp(i_E)
            disp(i_imm_freq)
            
            save(sprintf('%s',name_file))
            toc
        end
    end
    timeNeeded1=toc;
    disp(timeNeeded1);
else
    load(sprintf('%s',name_file))   
end
