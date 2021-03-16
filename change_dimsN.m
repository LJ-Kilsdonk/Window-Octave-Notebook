function change_dimsN(n_phenotypes,n_lineages)
%CHANGE_DIMSN changes the dimensions of N in the model
% changes the number of phenotypes a to n_phenotype and lineages n1 to
% n_lineage in the model. If either n_phenotype or n_lineage equals [], the
% respective dimension is not changed. 
%
% This function can affect the global parameters: MY_SETTINGS.Nelements,
% the dimensions of N, a, n1, and it also calls updateQandd to update the
% global parameters d and Q. 
%
% change_dimsN(n_phenotypes,n_lineages)
%
% See also: usemodel, makeQPheno, ode_ColonizationDynamics_standard,
% ode_ColonizationDynamics_stochastic_finite

%
global a n1 N MY_SETTINGS

% making sure that the inital "a" and n1 equal the number of rows and
% columns in N.
a = size(N,1);
n1 = size(N,2);

% shrink N or grow N with extra columns/rows filled with zeros
if ~isempty(n_phenotypes)
    if n_phenotypes > a
        N = [N; zeros(n_phenotypes - a,size(N,2))];
    else
        N = N(1:n_phenotypes,:);
    end
    a = n_phenotypes;
end
if ~isempty(n_lineages)
    if n_lineages > n1
        N = [N zeros(size(N,1), n_lineages - n1)];
    else
        N = N(:,1:n_lineages);
    end
    n1 = n_lineages;
end

% update Q and d, as they depend on the dimensions of N.
updateQandd;
% update number of elements in the model
MY_SETTINGS.Nelements = numel(N);
