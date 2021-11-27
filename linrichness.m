function linrichness(list_imm_freq, listE, cutoff,N_end1, N_end2)
%LINRICHNESS  Plot errorbar figures for lineage richness
% 
%    If N_end2 is provided, it is assumed to be without mutations, and both
%    the lineage richness in N_end1 and N_end2 are plotted together.
% 
%    LINRICHNESS(list_imm_freq, listE, cutoff,N_end1, N_end2) where
%    list_imm_freq is a list of the immigration frequencies used, listE is
%    the list of E values used, cutoff is the threshold in density above
%    which a lineage is included in the lineage richness count, and N_end1
%    and N_end2 are arrays of the final population density of each lineage
%    as given by [N_end,~]=loopmulti(...).
%    This function only works with one or two E values (i.e. length(listE)
%    < 3). 
% 
%   See also: loopmulti, linrichness
% 

if nargin > 4
    n_models = 2; %number of models used in plot
else
    n_models = 1;
end

for modeli=1:n_models
%    Model settings
    if modeli==1
        N_end = N_end1;
        Lstyle = {'-','-'};
        L_W = 2.8;
        colorsforplot={[196 33 98]/255;[57 151 135]/255};
        imm_freq = {list_imm_freq;list_imm_freq}; % only works for length(listE) < 3
    else
        N_end = N_end2;
        Lstyle = {'--','-.'};
        L_W = 1.5;
        colorsforplot={[0 0 0];[0 0 0]};
        imm_freq = {list_imm_freq + 0.0125;list_imm_freq - 0.0125}; % only works for length(listE) < 3
    end
%   calculate values for plots
    linRichA = squeeze(sum(N_end > cutoff, 1)); % lineage_richness(i_replicate, i_E, i_imm_freq)
    meanLR = squeeze(mean(linRichA,1)); % mean lineage richness(i_E, i_imm_freq)
    stdLR  = squeeze(std(linRichA,1)); % standard deviation lineage richness (i_E, i_imm_freq).
    
    for i_E= 1:length(listE)
        errorbar(imm_freq{i_E}, meanLR(i_E,:), stdLR(i_E,:));
%         errorbar(imm_freq{i_E}, meanLR(i_E,:), stdLR(i_E,:), ...
%             'LineStyle',Lstyle{i_E},'Color',colorsforplot{i_E},'LineWidth',L_W);      
        hold on
        xlabel('immigration frequency')
        ylabel('mean lineage richness')
    end
    set(gca,'FontSize',14);
    set(gca,'FontName','Ariel')
end

