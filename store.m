function store
%STORE  Sets the start values for state variable N to the final values of the
%   previous run, i.e. my_N(end,:)
%   
%   See also: runM, runMgrow

%   
global my_N N a
N = reshape(my_N(end,:), a, size(my_N,2)/(a));
