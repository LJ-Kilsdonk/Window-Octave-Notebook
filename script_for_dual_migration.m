%% Dual migration setup
% running these codes creates the results for the dual migration setup
% (both for the main text and the appendix). 
%
% For the simulations that determine the window of opportunity of
% opportunity for establishment and the simulation to determine adaptation
% rates, simulation output has already been saved. To not load the old
% simulation but redo them, change isredo to true in the section you want
% to redo.
%
% 
% 
% 
%% MULLER PLOTS
% MULLER PLOTS AS IN MAIN TEXT (examples of pre-equilibrium dynamics in runs)
%
% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = @ode_ColonizationDynamics_standard;
loadmodel

% set number of lineages to 2
n1 = 2;
change_dimsN([],n1);

% Run simulations with 2 immigration events and generate muller plots
%
% with invasion at t=5
figure;
dualmigrunandplot(15,14,5,300)

% with invasion at t=11
figure;
dualmigrunandplot(15,14,11,300)

% with invasion at t=20
figure;
dualmigrunandplot(15,14,20,300)

% 
% for saving use :print('NAME', '-dpng', '-r300') this gives high
% resolution files.


%% MULLER PLOTS VS POPULATION SIZE (WITH FINITE MODEL; SEE APPENDIX)
%
% 
% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = @ode_ColonizationDynamics_stochastic_finite;
loadmodel

% This random seed was chosen because it gives two contrasting examples in
% the two runs with K=10^7. 
rng(18)             

% set number of lineages to 2, and number of phenotypes to 10
n1 = 2;             % number of lineages
a  = 10;            % number of phenotypes
change_dimsN(a,n1); 

% colormap to use for the figures
MAP=jet(a);         % standard jet as color map
% MAP=gray(a+1);    %grayscale alternative
% MAP=MAP(2:end,:);

% Run simulations with 2 immigration events and generate muller plots
%
% run simulation with K=1e7:
fig_h1 = figure;
K=1e7;  
dualmigrunandplot(07,06,20,400,'colMAP',MAP,'whiteSpaceColor',[0 0 0])
title(sprintf('K=%g',K))

% run simulations again with K=1e7 (depending on random seed this can give
% very different results due to stochasticity)
fig_h2 = figure;
K=1e7;
dualmigrunandplot(07,06,20,400,'colMAP',MAP,'whiteSpaceColor',[0 0 0])
title(sprintf('K=%g',K))

% run simulation with K=1e6
fig_h3 = figure;
K=1e6;
dualmigrunandplot(07,06,20,500,'colMAP',MAP,'whiteSpaceColor',[0 0 0])
title(sprintf('K=%g',K))


% run simulation with K=1e15
fig_h4 = figure;
K=1e15;
dualmigrunandplot(07,06,20,250,'colMAP',MAP,'whiteSpaceColor',[0 0 0])
title(sprintf('K=%g',K))

% run simulation with K=1e4
fig_h5 = figure;
K=1e4;
dualmigrunandplot(07,06,20,1000,'colMAP',MAP,'whiteSpaceColor',[0 0 0])
title(sprintf('K=%g',K))

% Dual migration run with main text model, but with a=10 like the finite
% model; this is the infite population size version of the previous runs:
%
% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = @ode_ColonizationDynamics_standard;
loadmodel


% set number of lineages to 2, and number of phenotypes to 10
n1 = 2;             % number of lineages
a  = 10;            % number of phenotypes
change_dimsN(a,n1);

% run simulations with K = infinity
fig_h6 = figure;
dualmigrunandplot(07,06,20,250,'colMAP',MAP,'whiteSpaceColor',[0 0 0])
title('K=\infty');

%
%NOTE: for saving m?ller plots, use:
%print('NAME', '-dpng', '-r300') %this gives 300dpi resolution files.


%%
%----------WINDOW OF OPPORTUNITY PLOTS-------------------------------------
%% Standard main-text Window of opportunity for establishment
%

% settings to load data from previous simulations:
isredo=false;
% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel;

% set number of lineages to 2
n1 = 2;
change_dimsN([],n1);



% default behavior if the file to load from/save to doesn't exist, yet
% isredo is false, is to do the simulations anyway and save them under
% filename.
filename = 'WindowStandard.mat';% name of file to load from or save to
if isredo || ~(exist(filename,'file') == 2)
    % Find window of opportunity for all founder and invader phenotype
    % combos and store as matrix Wout:
    Wout = WindowSearch(0.05,500);
    save(filename,'Wout')
else
    load(filename,'Wout')
end

% plot the window of opportunity
plot_window(Wout);
 
%% window with delta death rate / mutations a = 10; dmin = 0; dmax = 18/19
%
%

% settings to load data from previous simulations:
isredo=false;
% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel;

% set number of lineages to 2 and set alternative values
n1 = 2;         % number of lineages
a = 10;         % number of phenotypes
d_min = 0;      % min. d(competition-independent death rate)
d_max = 18/19;  % max. d(competition-independent death rate)
d_diff = d_max - d_min; 
% change dimensions of N and update d and Q
change_dimsN(a,n1);


% default behavior if the file to load from/save to doesn't exist, yet
% isredo is false, is to do the simulations anyway and save them under
% filename.
filename = 'WindowAltDeathRatesv1.mat';% name of file to load from or save to
if isredo || ~(exist(filename,'file') == 2)
    Wout = WindowSearch(0.05,500);
    save(filename,'Wout')
else
    load(filename,'Wout')
end


% Plot the window of opportunity
%
% to allow for comparison with the standard window of opportunity (a=20),
% the same value for death rate gets the same color assigned
colorforlines=jet(20);
colorforlines=colorforlines(1:2:20,:);
plot_window(Wout,colorforlines);
xlim([0 125]) 

%% window with delta fitness / mutations a = 20; dmin = 15/38; dmax = 34/38
%

% settings to load data from previous simulations:
isredo=false;
% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel;


% set number of lineages to 2 and set alternative values
n1 = 2;         % number of lineages
d_min = 15/38;  % min. d(competition-independent death rate)
d_max = 34/38;  % max. d(competition-independent death rate)
d_diff = d_max - d_min;
% change dimensions of N and update d and Q
change_dimsN([],n1);

% default behavior if the file to load from/save to doesn't exist, yet
% isredo is false, is to do the simulations anyway and save them under
% filename.
filename = 'WindowAltDeathRatesv2.mat';
if isredo || ~(exist(filename,'file') == 2)
    Wout = WindowSearch(0.05,1000);
    save(filename,'Wout')
else
    load(filename,'Wout')
end

% Plot the window of opportunity
%
% to allow for comparison with the standard window of opportunity,
% the same value for death rate (in distance from optimum, should get the
% same color assigned.
colorforlines = jet(40);
% for only half the lines there is a line in the standard window of
% opportunity plot, with the same death rate (in distance from optimum
% (=minimum)). The rest is exactly in between; to allow for comparison,
% these are turned not plotted.
line_indeces = 1:2:19;
colorforlines = colorforlines(line_indeces,:);
plot_window(Wout,colorforlines,line_indeces);
xlim([0 500]) 


%% window with effect of propaguel size 10^-6
%
% settings to load data from previous simulations:
isredo=false;
% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel;

% set number of lineages to 2 and set alternative propagule size
n1 = 2;             % number of lineages
propsize = 1e-6;    % propagule size
change_dimsN([],n1);

% default behavior if the file to load from/save to doesn't exist, yet
% isredo is false, is to do the simulations anyway and save them under
% filename.
filename = 'WindowLowPropsize.mat';
if isredo || ~(exist(filename,'file') == 2)
    Wout = WindowSearch(0.05,750);
    save(filename,'Wout')
else
    load(filename,'Wout')
end

% Plot the window of opportunity
plot_window(Wout);
xlim([0 250])

%% window with effect of mutation rate 10^-4
%
% settings to load data from previous simulations:
isredo=false;
% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel;

% set number of lineages to 2 and set alternative mutation rate
n1 = 2;             % number of lineages
mut = 1e-4;         % mutation rate
change_dimsN([],n1);

% default behavior if the file to load from/save to doesn't exist, yet
% isredo is false, is to do the simulations anyway and save them under
% filename.
filename = 'WindowMutHigh.mat';
if isredo || ~(exist(filename,'file') == 2)
    Wout = WindowSearch(0.05,750);
    save(filename,'Wout')
else
    load(filename,'Wout')
end

% Plot the window of opportunity
plot_window(Wout);
xlim([0 250])

%% window with effect of mutation rate 10^-6
% settings to load data from previous simulations:
isredo=false;
% to redo the simulation; replace previous line with:
% isredo = true;

% load the standard (non-stochastic) model including all the default
% parameter settings  (and make them global in the command line)
modelname = 'ode_ColonizationDynamics_standard';
loadmodel;


% set number of lineages to 2 and set alternative mutation rate
n1 = 2;             % number of lineages
mut = 1e-6;         % mutation rate
change_dimsN([],n1);

% default behavior if the file to load from/save to doesn't exist, yet
% isredo is false, is to do the simulations anyway and save them under
% filename.
filename = 'WindowMutLow.mat';
if isredo || ~(exist(filename,'file') == 2)
    Wout = WindowSearch(0.05,750);
    save(filename,'Wout')
else
    load(filename,'Wout')
end

% Plot the window of opportunity
plot_window(Wout);
xlim([0 250])



%---------------ADAPTATION RATE--------------------------------------
%% EFFECT FINITE POPULATION SIZE ON ADAPTATION RATE
%
% settings to load data from previous simulations:
isredo=false;
% to redo the simulation; replace previous line with:
% isredo = true;

% determine adaptation rates and produce plot
adaptation_popsize(isredo);





