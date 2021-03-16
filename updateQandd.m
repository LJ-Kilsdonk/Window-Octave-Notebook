function updateQandd
%UPDATEQANDD initializes/updates Q and d
%
%   Initializes the mutation matrix Q and the resource-independent
%   death-rate d
% 
%   For simulations without mutations (i.e. mut=0) instead of letting Q be
%   the identity matrix, it is set to 1. which is faster.
%   If in addition to mut=0, a=1, the function assumes the faster vector-N
%   version is being used. In this case, the real number of possible
%   phenotypes has to be defined by the global variable a_orig. The
%   vector-N version uses the fact that without mutations each lineages
%   (i.e. column) only contains one non-zero phenotype (row). Instead of
%   the rows of N coding for phenotypes and the columns of N for lineages,
%   N is a column vector where each element corresponds to a lineage, and
%   that lineage has one phenotype (and thus one value for d). This d is
%   thus also a column vector. The phenotype of each lineage is the column
%   v_Ph, and has to be provided as a global variable for this version to
%   work.  
%   
%   needs to be run before simulation;
%   

global a mut d n1 E d_min d_diff Q v_Ph a_orig
% update/initialize the mutatation matrix Q:
if mut == 0
    Q = 1;
    
    if a == 1   % if true, use faster N-vector version
        % update/initialize the competition-independent death rate d
        v_j = abs(v_Ph(1:n1)' - E) / (a_orig-1); % match with environment per phenotype i
        d = v_j .* (d_diff) + d_min; % d of each lineage
        return
    end
else
    % mutants go up/down one index value in phenotype
    Q = diag(repmat(mut,a - 1,1),-1) + diag(repmat(mut,a - 1,1),1); 
    % add remaining non-mutating fractions
    Q = Q +   eye(a)-diag(sum(Q)); 
end

% update/initialize the competition-independent death rate d:
v_i = abs((1:a) - E) / (a-1); % match with environment per phenotype i
d_i = v_i .* (d_diff) + d_min; % d of each phenotype
d = repmat(d_i', [1 n1]); % matrix of d values with same dimensions as N
