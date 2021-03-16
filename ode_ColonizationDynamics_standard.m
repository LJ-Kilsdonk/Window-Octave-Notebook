function my_X2=ode_ColonizationDynamics_standard(t,my_X1)
%ODE_COLONIZATIONDYNAMICS_STANDARD  standard model as used in the main
%   text of Eco-evo dynamics of new populations. For a detailed description
%   of the model, please see [AmNat reference].
% 
%   [my_X2] = ode_ColonizationDynamics_standard(t,my_X1) returns the vector my_X2
%   with the different elements of the state variable N (a matrix). 
%   The input t (time) is not used, but is asked in order to conform with
%   solvers like ode45 and euler. my_X1 is a vector containing elements
%   of N at t, i.e. N(:). 
%      
%   The model uses global variables a d n1 Q r and MY_SETTINGS.Nelements
%   
global a d n1 Q r 
global MY_SETTINGS


% preallocate memory for output
my_X2 = zeros(MY_SETTINGS.Nelements, 1);
% convert input vector into a matrix
N = reshape(my_X1(1:MY_SETTINGS.Nelements), a, n1);
% Calculate the rates for the system of ordinary differential equations as
% and store it as a vector
my_X2(1:MY_SETTINGS.Nelements) = r.*(max(Q*N,0)) - r.*N.*repmat(sum(sum(N)),a,n1) - d.*N;