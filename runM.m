function runM(tspan, is_store_output) 
%RUNM   run model
%   runs the model and does not change the dimension of N.
%   Requires the model to be loaded to initialize random variables with the
%   function usemodel(modelname) or command-line code: loadmodel. 
% 
%   runM(tspan) with tspan = [t0 t1 ... t_final] runs the model and saves
%   tspan to my_t and the output of the model at those time steps to my_N.
%   my_N and my_t are global variables. my_N is a matrix with as row the
%   time, and as columns the elements of the state variable matrix N. 
%   
%   runM(tspan) with tspan = [t_final] uses the values 
%   [0:MY_SETTINGS.SolvOpt.MaxStep:t_final]. The default
%   MY_SETTINGS.SolvOpt.MaxStep = 0.1. 
%
%   runM(tspan) with tspan = [t0 t_final] uses as t values
%   [t0:MY_SETTINGS.SolvOpt.MaxStep:t_final]. 
%
%   runM(tspan, is_store_output) where if is_store_output is true, the last
%   values of the run (and row of my_N) are stored in N. The boolian
%   variable is_store_output is false by default.
%
%   See also runMgrow, euler, ode45, runMgrow, usemodel, loadmodel
%
global my_t my_N MY_SETTINGS t N
% set default
if nargin < 2
    is_store_output = 0;
end

% prevent error if dimensions of N have changed, but Q and d are not
% updated
updateQandd
% run the model
% [my_t, my_N] = feval(str2func(MY_SETTINGS.Solvname), str2func(MY_SETTINGS.modelname),...
%     tspan, N, MY_SETTINGS.SolvOpt);
 [my_t, my_N] = feval(MY_SETTINGS.Solver, MY_SETTINGS.model,...
     tspan, N, MY_SETTINGS.SolvOpt);
% set final N as new N
if is_store_output
    store
end
