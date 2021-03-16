function arrivalrank(N_end1, N_end2, list_imm_freq,listE)
%ARRIVALRANK    Plot figures for mean percentage of population with arrival
%   rank (j) comparing with and without mutation in one figure. Also
%   calculates the grand mean arrival rank with mutations and the effect
%   size (i.e. how much lower the grand mean arrival rank is with than
%   without mutations)
%   
%   ARRIVALRANK(N_end1, N_end2, list_imm_freq,listE) where N_end1 and
%   N_end2 are arrays of the final population density of each lineage as
%   given by [N_end,~]=loopmulti(...) with mutations and without mutations,
%   respectively. list_imm_freq is a list of immigraion frequencies used,
%   and listE is a list of E values used. 
% 
%   See also: loopmulti, linrichness
%   


% Test if dimensions are correct
if ~all(size(N_end1)==size(N_end2))
    error('the dimensions (n1max, m_repl, length(listE), length(list_imm_freq)) are not all identical for file_name1 and file_name2')
end
if size(N_end1,3)~=length(listE)
    error('length(listE) does not match the size(N_end,3); the number of E values tested')
end
if size(N_end1,4)~=length(list_imm_freq)
    error('length(list_imm_freq) does not match the size(N_end,4); the number of immigration frequencies tested')
end
% set dimensions
n1max  = size(N_end1,1);
n_E    = length(listE);
n_disp = length(list_imm_freq);

% Settings for plots
colorsforplot1 = mat2cell([[196 33 98]/255;[57 151 135]/255;hsv(max(n_E-2,0))],ones(max(n_E,2),1)); % Colors used in plots from N_end1
colorsforplot2 = mat2cell(repmat([150 150 150]/255,n_E,1),ones(n_E,1)); % Colors used in plots from N_end2

insetHpos = [0.18 0.62]; % List of horizontal positions of the insets
insetVpos = [0.75 0.45 0.15]; % List of vertical positions of the insets
insetWidth = 0.27; % Width of insets
insetHeight = 0.16;% Height of insets

% Calculating mean % of the population for different arrival ranks
% (averaging over the replicates)
% With mutations
pc_end1 = 100 * N_end1 ./ repmat(sum(N_end1,1),size(N_end1,1),1); % percentage of total density of each lineage within each replicate
PC_A1 = reshape(mean(pc_end1(:,:,:,:),2), [n1max,n_E,n_disp]); %mean % of the population over all replicates; dimensions: arrival_rank x E x imm_freq)

% Same without mutations
pc_end2 = 100 * N_end2 ./ repmat(sum(N_end2,1),size(N_end2,1),1);
PC_A2 = reshape(mean(pc_end2(:,:,:,:),2), [n1max,n_E,n_disp]);

% Max value for vertical axis limit (regular) plots
% rounded up to multiples of 5
PC_Amax = ceil(max([PC_A1(:); PC_A2(:)])/5)*5; 

for i_disp=1:n_disp
    % Maximum values for axis limits
    maxPc1         = max(max(PC_A1(:,:,i_disp))); % Max mean % of population with mutations irrespective of E
    maxPc2         = max(max(PC_A2(:,:,i_disp))); % Max mean % of population with mutations irrespective of E
    maxPc = ceil(max([maxPc1 maxPc2])/5)*5; % Max value for vertical axis limit inset
 
    for i_E = 1:n_E
        h1 = subplot(n_disp,n_E,i_E+n_E*(i_disp-1)); % Create subplot
        N_v1 = PC_A1(:,i_E,i_disp); % Vector of mean fractions in sims with mutations used in this iteration
        N_v2 = PC_A2(:,i_E,i_disp); % Vector of mean fractions in sims without mutations used in this iteration
        
        % The highest observed arrival rank
        maxArrivalRank1  = find(PC_A1(:,i_E,i_disp)>0,1,'last');
        maxArrivalRank2  = find(PC_A2(:,i_E,i_disp)>0,1,'last');
        
        % Indeces of used Ranks. These indeces are used to omit arrival
        % ranks that are higher than any observed arrival rank, instead 
        % of plotting them as a line at zero.
        Ranks1 = 1:maxArrivalRank1;
        Ranks2 = 1:maxArrivalRank2;        
        
        % Plot regular stair plots
        % To position tick in the middle of plateu, remove 0.5 (half the
        % size of each plateau (= horizontal part in the stair plot)
        stairs(Ranks1-0.5, N_v1(Ranks1), 'Color',colorsforplot1{i_E},'LineWidth',1.9)
        hold on
        stairs(Ranks2-0.5, N_v2(Ranks2), 'Color',colorsforplot2{i_E},'LineWidth',1.5)        
        xlim([0.5 size(N_v1,1)+0.5])
        ylim([0 PC_Amax])
        set(h1,'YTick',0:10:100);
        set(h1,'XTick',0:100:n1max);
        ylabel('mean % of population')
        xlabel('arrival rank (j)')
        set(gca,'FontSize',14);
        set(gca,'FontName','Ariel')
        
        % Plot insets ("zoomed-in versions") of stair plots
        fig1 = axes('Position',[insetHpos(i_E) insetVpos(i_disp) insetWidth insetHeight]);
        box on
        stairs(Ranks1-0.5, N_v1(Ranks1),'Color',colorsforplot1{i_E},'LineWidth',1.9);
        set(fig1,'XTick',0:10:60)
        set(fig1,'YTick',0:10:100)
        hold on
        stairs(Ranks2-0.5, N_v2(Ranks2),'Color',colorsforplot2{i_E},'LineWidth',1.5)
        xlim([0.5 50.5])
        ylim([0 maxPc])
        set(gca,'FontSize',14);
        set(gca,'FontName','Ariel')       
    end
end

% Grand mean arrival rank
% here we calculate the effect size (i.e. difference) in grand mean arrival
% rank
%
arrival_ranks = (1:n1max)';
% creates an array of same size as N_end1 (and N_end2) with the arrival ranks
arrival_ranksA = repmat(arrival_ranks,1,size(N_end1,2),n_E,n_disp);
% calculates frequency of each rank for every simulation
freq_by_Arrival_rank1 = (N_end1./repmat(sum(N_end1,1),n1max,1,1,1));
freq_by_Arrival_rank2 = (N_end2./repmat(sum(N_end2,1),n1max,1,1,1));
% The weighted mean of the arrival ranks (weighted by their frequency):
Wmean_Arrival_rank1 = sum(arrival_ranksA .* freq_by_Arrival_rank1, 1);
Wmean_Arrival_rank2 = sum(arrival_ranksA .* freq_by_Arrival_rank2, 1);
% Calculates the grand mean arrival rank (the mean arrival rank averaged
% over all replicates)
grand_mean_arrival_rank_1 = mean(Wmean_Arrival_rank1,2);
grand_mean_arrival_rank_2 = mean(Wmean_Arrival_rank2,2);
% The effectSize of the change in settings between with and without
% mutations is calculated using the difference of these grand means:
effectSizeD = squeeze(grand_mean_arrival_rank_2 - grand_mean_arrival_rank_1);
% display results
disp('effecSizeD:')
disp(effectSizeD)
disp('GrandMeanArrivalRank with mutations:')
disp(squeeze(grand_mean_arrival_rank_1))