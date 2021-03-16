function my_X2=ode_ColonizationDynamics_stochastic_finite(t,my_X1)
%ODE_COLONIZATIONDYNAMICS_STOCHASTIC_FINITE     stochastic model as used in
%   the appendix of Eco-evo dynamics of new populations. For a detailed
%   description of the model, please see [AmNat reference].
%
%   [my_X2] = ode_ColonizationDynamics_stochastic_finite(t,my_X1) returns
%   the vector my_X2 with the different elements of the state variable N (a
%   matrix).  The input t (time) is not used, but is asked in order to
%   conform with solvers like ode45 and euler. my_X1 is a vector containing
%   elements of N at t, i.e. N(:). 
%      
%   The model uses global variables a h K d n1 Q r and MY_SETTINGS.Nelements
%   make sure that MY_SETTINGS.SolvOpt.MaxStep == h
%   
% 
% 
global a h K d n1 Q r
global MY_SETTINGS
% preallocate memory for output
my_X2=zeros(MY_SETTINGS.Nelements,1);
% convert input vector into a matrix
N = reshape(my_X1(1:MY_SETTINGS.Nelements),a,n1);
% Calculate the rates for the system of ordinary (stochastic differential
% equations as and store it as a vector.
my_X2(1:MY_SETTINGS.Nelements) = ...
    (poissrnd(h.*(max(Q*N,0)./sum(sum(max(Q*N,0)))) .*(K.*sum(sum(r.*max(Q*N,0))))) ./(K.*h)) ...
    - r.*N.*repmat(sum(sum(N)),a,n1) - d.*N;

