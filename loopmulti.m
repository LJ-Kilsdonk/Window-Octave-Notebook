function [N_end,MultiM_opt] = loopmulti(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file,MultiM_opt)
%LOOPMULTI  Main loop for runs with multiple migrations
%   
%   Loops over different settings and replicates, each iteration it calls
%   multi_migration_run to perform a simulation with multiple immigration
%   events, for which it stores the final population density of each
%   lineage. 
% 
%   [N_end,MultiM_opt] = LOOPMULTI(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file,MultiM_opt)
%   where N_end is the final population density of each lineage with
%   dimensions max_number_of_lineages x m_rep x length(listE) x
%   length(list_i_imm_freq), MultiM_opt is a structure containing the
%   settings used for the simulations and can be changed to not use default
%   settings, m_repl is the number of replicates, timeLrun the duration of
%   a run, list_imm_freq a list of the immigration frequencies to be used,
%   listE a list of the E values to be tested, isredo a boolean operator
%   that states if the simulations are redone, even if name_file already
%   exists. If false and the file exists, the output data is loaded from
%   the file.
%   Default maximum number of lineage is 350. To change this, change
%   MultiM_opt.n1max  
% 
%   See also: preequilloopmulti, multi_migration_run

global propsize E a mut
% Default behavior if the file to load from/save to doesn't exist, yet
% isredo is false, is to do the simulations anyway and save them under
% filename.
if isredo || ~(exist(name_file,'file') == 2)
    mut = getknownfield(MultiM_opt,'mut',mut);
    MultiM_opt.warnings =[]; % not used in current version; but initializes MultiM_opt if not provided.
    
    % create probability distribution of phenotype frequencies among migrants
    amin1 = a - 1;
    S = (factorial(amin1)./(factorial(0:amin1).*factorial(amin1- (0:amin1))))./(2^amin1); %phenotype frequencies among migrants
    MultiM_opt.S = S;
    MultiM_opt.propsize = getknownfield(MultiM_opt,'propsize',propsize);
    
    % Previous versions allowed for grouping several immigrants under the
    % same 'lineage group'. Here each migrant corresponds to one lineage.
    % Therefore the column index of each lineage is simply the arrival rank
    % of that lineage.
    MultiM_opt.n1max = getknownfield(MultiM_opt,'n1max',350);
    MultiM_opt.lin_j=(1:MultiM_opt.n1max)';
    
    % minimum stepsize between two invasions.
    MultiM_opt.inv_stepsize = getknownfield(MultiM_opt,'inv_stepsize',0.1);
    inv_stepsize = MultiM_opt.inv_stepsize;
    
    % Allocate memory for loop
    N_end=zeros(MultiM_opt.n1max,m_repl,length(listE),length(list_imm_freq));
    
    % Main Loop
    tic

    disp('paramters that will be tested in loopmulti:')
    disp('E:')
    disp(listE)
    disp('immigration frequency:')
    disp(list_imm_freq)
    disp('number of repetitions:')
    disp(m_repl)
    for i_E = 1:length(listE)
        E = listE(i_E);
        for i_imm_freq = 1:length(list_imm_freq)
            imm_freq = list_imm_freq(i_imm_freq);
            
            % display settings to monitor progress
            disp('start testing parameter values:')
            disp('E=')
            disp(E)
            disp('immigration frequency=')
            disp(imm_freq)            
            disp('printing repetitions:')
            for i_rep = 1:m_repl
                % Store information unique to current loop (for debugging)
                MultiM_opt.indIter = [E imm_freq i_rep];
                % Run model and store final population density of each
                % lineage.
                [N_end(:,i_rep,i_E,i_imm_freq),] = multi_migration_run(S,imm_freq,timeLrun,inv_stepsize,MultiM_opt);
                % display repetition to monitor progress
                disp(i_rep);
            end
            toc            
            save(sprintf('%s',name_file))
        end
    end
else
    load(sprintf('%s',name_file))
end

