function enterfounder(phenotype,propagule_density)
%ENTERFOUNDER  Sets the state variable such that it contains that only
%   the specified phenotype.
%
%   ENTERFOUNDER(phenotype, propagule_density) where possible values for
%   phenotype are in the set 1:a. Default phenotype is the least adapted
%   phenotype (phenotype = a). Default propagule_density = propsize, and if
%   propsize is undefined, it equals 0.001.
%
%

%
global N propsize a n1
% if no propagule_density is provided, use the global variable propsize
if nargin<2
    if ~isempty(propsize)
        propagule_density = propsize;
        disp('propagule_size=')
        disp(propagule_density)
    else
        propagule_density = 0.001;
        disp('no propagule_size defined, propagule_size set to 0.001')
    end
end

if nargin>0 && phenotype>a
    warning('phenotype was set at %d, but ther are only %d phenotypes, phenotype is set to %d',phenotype,a,a);
    phenotype=a;
end

N=zeros(a,n1);
if nargin>0
    N(phenotype,1) = propagule_density;
end
