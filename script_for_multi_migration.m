%% script for multi-migration setup 
% run this script to create the results for the multiple migration setup
% (both for the main text and the appendix). 
%
% Simulation output has already been saved. To not load the old %
% simulation but redo them, change isredo to true in the section you want
% to redo. 
%

%
%% Equilibrium lineage richness and arrival rank (as used in main text):
%
%
% 1. WITH evolution of lineages (mut>0).
clear

% settings to load data from previous simulations:
isredo = false;
% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel;


% simulation settings:
rng(23);                     % random seed used for simulations in paper
m_repl = 200;                 % number of replicates; decrease m_repl to reduce simulation time.
listE = [1 11];              % list of E values tested
timeLrun=500;               % time length of the simulations (in generations)
list_imm_freq=[0.05 0.2 0.8];    % list of immigration frequency
name_file = 'multi_migration_setup_evolLin.mat'; % name of file to load from or save to
MultiM_opt.mut = mut;       % the structure multiM_opt is used to pass on settings to the main loop
% To determine max number of invasions, the
% following test can be run. The outcome is 500:
% n1max=n1maxtest(list_imm_freq(end),timeLrun,0.1, m_repl); % calculate maximum expected n1
% multiM_opt.n1max = ceil(n1max/50)*50; % round upwards to a nice round number (multiples of 50)
MultiM_opt.n1max = 500; % maximum number of invasions.

% Plot settings
% Percentage of propagule density used as cutoff population density for
% inclusion in lineage richness count. 
cutoff_as_perc_of_propsize = 105;

% run simulations or load data from previous simulations
[N_end_mut,MultiM_opt_mut] = loopmulti(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file,MultiM_opt);


% 2. WITHOUT evolution of lineages (mut=0).
% settings to load data from previous simulations:
isredo=false;
% to redo the simulation; replace previous line with:
% isredo = true;

% reload the standard (non-stochastic) model including all the default
% parameter settings (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel

% changes in simulation settings:
rng(33);                    % random seed used for simulations in paper
mut = 0;                    % no mutations
MultiM_opt.mut = mut;       % the structure multiM_opt is used to pass on settings to the main loop
name_file = 'multi_migration_setup_noMut.mat'; % name of the file to save / load from the simulations without mutations
% run simulations or load data from previous simulations
[N_end_nomut, MultiM_opt_nomut] = loopmulti(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file,MultiM_opt);

% plot results for BOTH WITH and WITHOUT evolution of lineages:
f_h = figure;
arrivalrank(N_end_mut, N_end_nomut, list_imm_freq,listE);
% lineage richness needs a cutoff value for the population density,
% otherwise all 
if MultiM_opt_mut.propsize == MultiM_opt_nomut.propsize
    propdensity = MultiM_opt_mut.propsize; % propagule density
    cutoff = propdensity .* cutoff_as_perc_of_propsize/100;
else
    error('propagule density differs between simulations with vs. without mutations')
end
f_h2=figure;
linrichness(list_imm_freq, listE, cutoff, N_end_mut, N_end_nomut);



%% Pre-equilibrium simulations and figures as used for main text
% 

clear
% settings to load data from previous simulations:
isredo=false;

% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel;

% simulation settings:
rng(42);                % random seed used for simulations in paper
timeLrun=750;           % time length of the simulations (in generations)
m_repl=200;             % number of replicates; decrease m_repl to reduce simulation time.
list_imm_freq=[0.05 0.15];   % immigration frequency
listE= [1 11];          % list of E values tested
name_file = 'PreEquil_standard.mat'; % name of the file to save simulations to or load them from

% run simulations or load data from previous simulations (here we do both
% with and without mutations in the same loop).
[linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,MultiM_opt] = ...
    preequilloopmulti(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file);

% plot results
preEquilFigs(linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,list_imm_freq,listE,MultiM_opt);




%% Pre-equilibrium runs with main model for appendix 
% with infinite population size, USED FOR COMPARISON IN THE APPENDIX. It
% has only 10 instead of 20 phenotypes like the results in the main text.
%
clear
% settings to load data from previous simulations:
isredo=false;

% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel;

% simulation settings:
rng(44)                 % random seed used for simulations in paper
timeLrun=750;           % time length of the simulations (in generations)
list_imm_freq=[0.05 0.15];   % immigration frequency
listE= [1 6];           % list of E values tested
m_repl=200;             % number of replicates; decrease m_repl to reduce simulation time.
name_file = 'PreEquil_standard_a10.mat'; % name of the file to save simulations to or load them from

% change model settings for using 10 phenotypes
a = 10;
change_dimsN(a,[]);

% run simulations or load data from previous simulations (here we do both
% with and without mutations in the same loop).
[linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,MultiM_opt] = ...
    preequilloopmulti(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file);
% plot results:
preEquilFigs(linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,list_imm_freq,listE,MultiM_opt);


%% pre-equilibrium run with finite model used for the appendix K = 10^4
% uses the (stochastic) finite population size model. Includes with and
% without evolution
clear
% settings to load data from previous simulations:
isredo=false;

% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings (and make them global in the command line):
modelname = 'ode_ColonizationDynamics_stochastic_finite';
loadmodel;

% simulation settings:
rng(43)                 % random seed used for simulations in paper
timeLrun=750;           % time length of the simulations (in generations)
m_repl=200;             % number of replicates; decrease m_repl to reduce simulation time.
list_imm_freq=[0.05 0.15];   % immigration frequency
listE= [1 6];           % list of E values tested
name_file = 'PreEquil_stochastic_finite_a10.mat'; % name of the file to save simulations to or load them from

% run simulations or load data from previous simulations (here we do both
% with and without mutations in the same loop).
[linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,MultiM_opt] = ...
    preequilloopmulti(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file);
% plot results:
preEquilFigs(linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,list_imm_freq,listE,MultiM_opt);

%% pre-equilibrium run with finite model used for the appendix K = 2*10^4
% uses the (stochastic) finite population size model. Includes with and
% without evolution
clear
% settings to load data from previous simulations:
isredo=false;

% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings (and make them global in the command-line):
modelname = 'ode_ColonizationDynamics_stochastic_finite';
loadmodel;

% simulation settings:
rng(45)                 % random seed used for simulations in paper
timeLrun=750;           % time length of the simulations (in generations)
m_repl=200;             % number of replicates; decrease m_repl to reduce simulation time.
list_imm_freq=[0.05 0.15];   % immigration frequency
listE= [1 6];           % list of E values tested
K = 20000;              % population size
name_file = 'PreEquil_stochastic_finite_a10_K20000.mat'; % name of the file to save simulations to or load them from

% runs simulations or load data from previous simulations (here we do both
% with and without mutations in the same loop).
[linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,MultiM_opt] = ...
    preequilloopmulti(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file);
% plot results:
preEquilFigs(linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,list_imm_freq,listE,MultiM_opt);


%% pre-equilibrium run with finite model used for the appendix K = 10^5
% uses the (stochastic) finite population size model. Includes with and
% without evolution
clear
% settings to load data from previous simulations:
isredo=false;

% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings (and make them global in the command-line):
modelname = 'ode_ColonizationDynamics_stochastic_finite';
loadmodel;
% simulation settings:
rng(46)                 % random seed used for simulations in paper
timeLrun=750;           % time length of the simulations (in generations)
m_repl=200;             % number of replicates; decrease m_repl to reduce simulation time.
list_imm_freq=[0.05 0.15];   % immigration frequency
listE= [1 6];           % list of E values tested
K = 1e5;                % population size
name_file = 'PreEquil_stochastic_finite_a10_K100000.mat'; % name of the file to save simulations to or load them from

% runs simulations or load data from previous simulations (here we do both
% with and without mutations in the same loop).
[linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,MultiM_opt] = ...
    preequilloopmulti(m_repl,timeLrun,list_imm_freq,listE,isredo,name_file);
% plot results:
preEquilFigs(linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,list_imm_freq,listE,MultiM_opt);
