function v = getknownfield(s, f, d)
%GETKNOWNFIELD  Get field f from struct s, or else yield default d.
%
%   See also: odeget

%   Script is copied from the local function in the Matlab build-in odeget
%   function file.

if isfield(s,f)   % s could be empty.
    v = s.(f);
    if isempty(v)
        v = d;
    end
else
    v = d;
end