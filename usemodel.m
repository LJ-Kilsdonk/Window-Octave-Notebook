function usemodel(modelname)
%USEMODEL Loads the model
% This function creates global variables and initializes the specified
% model by giving these global variables the values belonging to that
% model.
%
% all parameters in the model (and propsize) are stored as a global
% variable. Other model specific settings, like the solver to be used and
% the model name, are stored in the structure: MY_SETTINGS
% 
% the two possible models to use are ode_ColonizationDynamics_standard and
% ode_ColonizationDynamics_stochastic_finite. The argument modelname can 1)
% contain the whole model name as a string, the model name as a string with
% the "ode_" or the "ode_ColonizationDynamics_" part removed, or it can be
% a function handle to the model.
%
% usemodel(modelname) 

%
global MY_SETTINGS propsize
global r mut n1 a E d d_min d_diff Q K h N

%handling different types of input for modelname. A function handle is
%allowed, but also a string containing the either the whole name of the
%function, the name without "ode_" at the start, or the name without
%"ode_ColonizationDynamics_" at the start.
if isequal(class(modelname), 'function_handle')
    MY_SETTINGS.modelname = func2str(modelname);
elseif ischar(modelname)
    if isequal(modelname(1:25),'ode_ColonizationDynamics_')
            MY_SETTINGS.modelname = modelname;
    elseif isequal(modelname(1:21),'ColonizationDynamics_')
        MY_SETTINGS.modelname=['ode_',modelname];
    else        
        MY_SETTINGS.modelname=['ode_ColonizationDynamics_',modelname];
    end    
else
    error('modelname has class %s. Only funcation_handle or string allowed',class(modelname))
end
MY_SETTINGS.model = str2func(MY_SETTINGS.modelname);


%---------default MODEL SPECIFIC parameter values/settings----------------
switch MY_SETTINGS.modelname
    case 'ode_ColonizationDynamics_standard'
        % Setting default parameter values specific to model
        a = 20;            % # of phenotypes
        
        % Default settings solver, specific to model
        MY_SETTINGS.Solvname = 'ode45';     %solver to be used
        MY_SETTINGS.Solver  = str2func(MY_SETTINGS.Solvname);
        MY_SETTINGS.SolvOpt.AbsTol = 1e-08;     
        MY_SETTINGS.SolvOpt.RelTol = 5e-06;
    case 'ode_ColonizationDynamics_stochastic_finite'
        % Setting default parameter values specific to model
        a = 10;            % # of phenotypes
        h = 0.1;           % size of time steps.
        K = 10000;         % population size
                
        % Default settings solver, specific to model
        MY_SETTINGS.Solvname = 'euler'; %solver to be used
        MY_SETTINGS.Solver  = str2func(MY_SETTINGS.Solvname);
        MY_SETTINGS.SolvOpt.MaxStep = h;
    otherwise
        error('unkown function (name): %s',MY_SETTINGS.modelname);
end


%---non_model_specific default/initial parameter/settings/variable values--
n1 = 100;          % #of lineages;
E = 1;             % best phenotype for environment in target patch
r = 1;             % maximum growth rate
mut = 1E-5;        % mutation rate
propsize = 0.001;  % propagule size
d_min = 0;         % min. d(competition-independent death rate)
d_max = 1;         % max. d(competition-independent death rate)
d_diff = d_max - d_min;
% d: competition-independent death rate; set with updateQandd
% Q: the mutation matrix; set with the function updateQandd
% Q and d are based on other parameters, and before each run
% they are recalculated using the function updateQandd.
updateQandd;

% The initial value state variable.
% Matrix of population densities with rows indexing phenotype and
% columns indexing (migrant) lineages):
N = [rand([a 1]).*0.01 zeros(a,n1-1)];

% Default settings of model
MY_SETTINGS.Nelements = numel(N);
