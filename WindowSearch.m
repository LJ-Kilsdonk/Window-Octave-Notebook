function W = WindowSearch(T,time_lengthrun)
%WINDOWSEARCH  Finds the window of opportunity for establishment
%   Searches for the window of opportunity for establishment using the
%   build-in matlab function fzero. 
%
%   W = WindowSearch(T,time_lengthrun)
%   where W is the lag in arrival time for all combination of founder and
%   invader phenotype, T is the frequency of the invader used as threshold
%   for success, and time_lengthrun is the time length of each run in
%   generations. 
%
global a

T_f = 1 - T; % minimal frequency founding lineage
optfsearch = optimset('Display','iter');

% memory allocation for loop
W = zeros(a,a);
Rlimit = time_lengthrun * 0.5; % the (initial tried) upper limit for the search function

%MAIN LOOP
% steps in loop:
% 1. is establishment of invader possible at all? (assumption: if it cannot
% establish when arriving simultaniously with the founder, establishment is
% impossible)
% 2. is the Rlimit the correct maximum window of opportunity (if not raise
% it)
% 3. find the lag in arrival time where T_f == frequency of the founder in
% equilibrium.
%
for i_PhC = 1:a
    for i_PhI = 1:a 
        frqFlin_sim = dualImm(0,i_PhC,i_PhI,time_lengthrun); %freqency of founding lineage with simultanious migration.
        if frqFlin_sim > T_f
            % Establishment never possible. Therefore we denote the window
            % of opportunity for establishment as 0. Note this does not
            % mean, invasion is possible when the lag in arrival time of the
            % invader equals 0. By definition the invader (i.e the 2nd
            % immigrant) has a positive non-zero lag in arrival time with
            % the founder (= colonizer = 1st-migrant). 
            W(i_PhC,i_PhI) = 0;
        else
            % if the Rlimit is too small it will be increased here
            while dualImm(Rlimit,i_PhC,i_PhI,time_lengthrun) < T_f
                Rlimit = Rlimit + time_lengthrun * 0.1;
                if Rlimit > time_lengthrun
                    disp(fprintf('phenotype colonist: %s',i_PhC))
                    disp(fprintf('phenotype invader: %s',i_PhI))
                    error('The window of opportunity is longer than 90% of time_length')
                end
            end
            % find the window of opportunity using fzero
            W(i_PhC,i_PhI) = fzero(@(WD) diff2T(WD,i_PhC,i_PhI,T_f,time_lengthrun), [0 Rlimit], optfsearch);
        end
        disp(i_PhC)
        disp(i_PhI)
    end
end

%----------------------------ToMinimize--------------------------
function ToMinimize = diff2T(WD,colonizer_pheno,invader_pheno,T_f,time_lengthrun)
%TOMINIMIZE     Function that is minimized.
%   Function whose absolute value needs to be minimized. This value is the
%   the frequency of the founder lineage minus the threshold for succesfull
%   establishment T_f.
%   ToMinimize = diff2T(WD,colonizer_pheno,invader_pheno,T_f,time_lengthrun)
frqFounderlineage = dualImm(WD,colonizer_pheno,invader_pheno,time_lengthrun);
ToMinimize = frqFounderlineage - T_f;


%-----------------Dual immigration run-------------------------
function frqFounderlineage = dualImm(WD,colonizer_pheno,invader_pheno,time_lengthrun)
%FRQFOUNDERLINEAGE  Dual immigration that returns the frequency of the
%   founder lineage
% 
%   frqFounderlineage = dualImm(WD,colonizer_pheno,invader_pheno,time_lengthrun)
%   function that performs a run with two migrants of given phenotypes with a
%   given lag in arrival time WD, phenotype of the first migrant
%   colonizer_pheno, phenotype of the 2nd migrant (invader_pheno), and length
%   of the run (time_lengthrun); and outputs the frequency in the population
%   of the descendants from the 1st migrant(founder): frqFounderlineage.
% 

global propsize N my_N
% Run model
enterfounder(colonizer_pheno,propsize)
% Run model upto the arrival of the 2nd immigrant (if WD=0, no run
% required)
if WD > 0
    runM([0 WD],1);
end
% Introduce 2nd immigrant
N(invader_pheno,2) = propsize;
% Continue run till end (i.e. time_lengthrun) and request also output one
% time step before end, which is used to verify that stablity has been
% reached.
runM([WD time_lengthrun-1 time_lengthrun],1); 
frqFounderlineage = sum(N(:,1),1)/sum(N(:));

% Test if N has converged on a stable equilibrium (note, the system
% will never reach complete stability, it's a limit. Therefore, 1e-5 is used
% as a tolerance level. If the change in population density is less than
% 1e-5, the system is practically in equilibrium
if sum(my_N(end,:),2) - sum(my_N(end-1,:),2) > 1e-5
    warning('simulations ended before stable equilibrium')
end