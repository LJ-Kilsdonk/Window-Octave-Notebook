function runMgrow(tspan) 
%runMgrow   Run model for the multi-migration setup.
%   Runs the model and aftewards adds one column to N, to allow for the
%   immigration of one new lineage.
%   Requires the model to be loaded to initialize random variables with the
%   function: "usemodel.m" or command-line code: "loadmodel.m". 
%   
%   runMgrow(tspan) with tspan = [t0 t1 ... t_final] runs the model and saves
%   tspan to my_t and the output of the model at those time steps to my_N.
%   my_N and my_t are global variables. my_N is a matrix with as row the
%   time, and as columns the elements of the state variable matrix N. The
%   last values of the run (and row of my_N) are stored in N.
%   
%   runMgrow(tspan) with tspan = [t_final] uses the values 
%   [0:MY_SETTINGS.SolvOpt.MaxStep:t_final]. The default
%   MY_SETTINGS.SolvOpt.MaxStep = 0.1. 
%
%   runMgrow(tspan) with tspan = [t0 t_final] uses as t values
%   [t0:MY_SETTINGS.SolvOpt.MaxStep:t_final]. 
%
%   
%
%   See also: runM, euler, ode45, usemodel, loadmodel

%   
global my_t my_N MY_SETTINGS t N n1

% run model
runM(tspan, true)

% add new column lineage
n1 = n1+1;
N(:,n1) = 0;
MY_SETTINGS.Nelements = numel(N);