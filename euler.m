function [tout,yout] = euler(odefun, tspan, y0, options)
%EULER  Solve system of ordinary diferential equations using the
%   Euler method with fixed time steps.
%
%   [tout,yout]=euler(odefun,tspan,y0,options) with tspan = [t0 t1 ...
%   tfinal] integrates the system of differential equations y' = f(t,y)
%   from time t0 to tfinal with initial conditions y0. tspan has to be
%   monotonic increasing or decreasing. odefun is a function handle. For
%   consistency with other solvers, like ode45, the odefun must request 2
%   input arguments: time (t) and state variable (y). options.maxstep is
%   the integration step size; by default it is 0.1. 
%
%   [tout,yout]=euler(odefun,tspan,y0,options) with tspan = tfinal
%   returns all calculated time steps: [0:options.Maxstep:tfinal].
%
%   [tout,yout]=euler(odefun,tspan,y0,options) with tspan = [t0 tfinal]
%   returns all calculated time steps: [t0:options.Maxstep:tfinal].
%
%   note: Because the euler function uses constant time-steps (=
%   options.Maxstep), it's advisable to use time steps in tspan that are a
%   multiples of options.Maxstep. Otherwise, linear interpolation is used
%   to approximate the value at a time points between two time points
%   calculated with the euler method.
%
%   See also: ode45

%   

% arguments for solver
delta = odeget(options,'MaxStep',0.1); % time steps in euler function, default = 0.1
tout = tspan(:); % time points used for output

% if only one value for tout, output gets all calculated time steps from 0
% to tout; If two values for tout, output gets all calculated time steps
% from t0 to tout.
if length(tout) <= 2
    if length(tout) == 1
        t0 = 0;
        tfinal = tout;
        tdir = 1;   % time direction: -1 = backward; 1 = forward
    elseif length(tout) == 2
        t0 = tout(1);
        tfinal = tout(end);  
        tdir = sign(tfinal - t0); % time direction
        delta = delta * tdir;
    end
    tout = t0:delta:tfinal;
    % to ensure that tfinal itself is in the output time steps
    if tdir * (tfinal - tout(end)) > 0
        tout(end + 1) = tfinal;
    end
else % i.e. if length(tout) >= 3
    t0 = tout(1);
    tfinal = tout(end);
    % allow backwards runs, and make sure tspan is monotonic (increasing or
    % decreasing):
    tdir = sign(tfinal - t0); % time directions: -1 = backward; 1 = forward
    if any( tdir*diff(tout) <= 0 )
        error('TspanNotMonotonic');
    end
    delta = delta * tdir;    
end
ntspan = length(tout); % number of time points for output


% time steps used internally
t_intern= t0 : delta : tfinal; 
% ensure the last time point in the output is included or at least inside
% the range of time values used internally.
if tdir * (tout(end) - t_intern(end)) > 0   % time direction changes sign
    t_intern(end + 1) = t_intern(end) + delta;
end
nt_intern_span = length(t_intern); % number of time point for internal use

% allocate memory for output
yout = zeros(ntspan, numel(y0)); % state variable values for output
yintern = zeros(nt_intern_span, numel(y0)); % state variable values used internally
% Note: output at every calculated timestep is saved internally. If ntspan
% is much smaller than nt_intern_span, a more efficient algorithm might be used. 


% THE MAIN LOOP
y = y0(:); % state variable
t = t_intern(1); % time used in model
yintern(1,:) = y.';
for i_tnew = 2: nt_intern_span 
    % use (system of) ode to find derivative y with respect to time
    f = feval(odefun, t, y); 
    % Euler method to find value y at next time step
    y = y + f .* delta;
    % update t to next time step
    t=t_intern(i_tnew);
    % save value y
    yintern(i_tnew, :) = y.';    
end

% if t_intern(i) = tout(i), interp1 gives yout(i,:) = yintern(i,:), else yout
% will be a linear interpolation of the two surrounding time points in
% tintern. e.g. if tout(4) = 0.15, t_intern(2)=0.1, and t_intern(3)=0.2,
% then yout(4) = yintern(2) + ( (yintern(3)-yintern(2)) / (0.2 - 0.1) ) *(0.15 - 0.1)
yout(1,:) = yintern(1, :); % as t0 = tout(1), t_intern(1) = tout(1), always.
yout(2:end,:) = interp1(t_intern, yintern, tout(2 : end));