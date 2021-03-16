function adaptation_popsize(isredo)
%ADAPTATION_POPSIZE   Adaptation time with different population sizes.
%   Tests the effect of population size on adaptation time with a = 10. The
%   function produces a boxplot based on 100 replicates, for each tested
%   population size: K = [1e4 2e4 logspace(5,15,11)], plus a horizontal
%   line to indicate the (non-stochastic) adaptation time of a run with
%   infinite population size. 

%   ADAPTATION_POPSIZE(isredo) if isredo is false, a previously saved file
%   with results is used.

global propsize K
if isredo || ~(exist('adaptation_popsize.mat','file') ==2)
    % Load model
    usemodel(@ode_ColonizationDynamics_stochastic_finite)
    % Change dimensions of N to 2 lineages.
    n1 = 2;
    propsize = 0.001;
    change_dimsN([],n1);
    
    % Set parameters
    K_list=[1e4 2e4 logspace(5,15,11)];
    m_rep=100;
    
    % Allocate memory for loop
    tim2totadpt=zeros(length(K_list),m_rep);
    
    % Main loop (with finite model).
    for i=length(K_list):-1:1
        K=K_list(i);
        disp(K);
        tim2totadpt(i,:) = findadaptationtime(m_rep);
    end
    
    % Same idea but now with the standard (infinite pop. size) model
    usemodel(@ode_ColonizationDynamics_standard)
    n1 = 2;
    a = 10;
    propsize = 0.001;
    change_dimsN(a,n1);
    tim2totadptINF = findadaptationtime(1);
    
    save('adaptation_popsize.mat')
else
    load('adaptation_popsize.mat')
end


% Plot results
figure;
horizontalBoxPosition = [4 log10(2*10^4) 5:15] - 3;
labelsBoxes = {'10^4','2x10^4','10^5','10^6','10^7','10^8','10^9','10^10',...
    '10^11','10^12','10^13','10^14','10^15'};
% To make exponents superscript the figure has to be edited further
% outside matlab. 
boxplot(tim2totadpt','labels',labelsBoxes,'Widths',0.1,'Positions',horizontalBoxPosition)
ylabel('time untill fully adapted')
xlabel('individuals at K')
hold on
line([0 14], [tim2totadptINF tim2totadptINF])
set(gca,'FontSize',14);
set(gca,'FontName','Ariel')


%-------------------------------------------------------------
function adapttime = findadaptationtime(m_rep)
%ADAPTTIME   Finds the time needed to adapt.
% The maximum density that can be reached in the model is 1. As there is 1)
% stochasticity, 2) even in the deterministic model this value is only
% reached in the limit (at t= infinity), and 3) due to mutations the
% population will never be 100% the most adapted phenotype, we use 0.995 as
% threshold population density for the most adapted phenotype, to call a
% population 'fully adapted'.
% 
global a propsize my_t my_N
adapttime = zeros(m_rep,1);
timepast=0;
for j=1:m_rep
    % Initialize state variable N
    enterfounder(a-1,propsize)
    
    % A priori it is not clear how many time-steps are minimally needed
    % for the population to adapt. To save time, it is best not to
    % strongly overestimate the adaptation time.
    % As estimation we use 120% of the adaptation time used in the
    % previous replicate. For the first replicate with a certain K
    % value, 500 time steps are used.
    % if the adaptation time is not sufficient 1000 time steps are
    % added.
    if timepast==0
        runM([0 500],1);
    else
        runM([0 timepast*1.2],1);
    end
    timepast=0;
    
    % To prevent getting stuck in a loop, after maximum 100 tries, we give up on finding the adaptationtime.
    for ij=1:100
        if any(my_N(:,1)>0.995)
            adapttime(j)= timepast + my_t(find(my_N(:,1)>0.995,1));
            timepast=adapttime(j);
            disp(timepast);
            break
        elseif sum(my_N(end,:))<0.001
            adapttime(j)=NaN;
            timepast=0;
            disp(timepast);
            disp('population went extinct')
            break
        end
        timepast=timepast + my_t(end);
        if ij==100
            error('no max adaptation reached');
        end
        runM([0 1000],1);
    end
end

