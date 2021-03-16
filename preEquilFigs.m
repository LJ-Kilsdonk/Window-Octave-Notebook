function preEquilFigs(linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,list_imm_freq,listE,MultiM_opt)
%PREEQUILFIGS   plot figures of mean lineage richness vs time and mean
%   arrival rank vs time. This function is used to visualize the output of
%   the pre-equilibrium simulations. 
%   
%   preEquilFigs(linR_Array,gmeanARank_Array, linR_Array_nomut,gmeanARank_Array_nomut,list_imm_freq,listE,MultiM_opt)
%   plots the output of preequilloopmulti. For input argument description
%   see preequilloopmulti.
% 
%   See also: preequilloopmulti
% 

% Settings for plots
colorsforplot1={[196 33 98]/255;[57 151 135]/255};
colorsforplot2={[0 0 0];[0 0 0]};
Lstyle1 = {'--';'-'};
Lstyle2 = {'--','-'};

% Determine axes limits for plots
x_up_limit = MultiM_opt.Roundedtime(end);
maxLinRich = max([linR_Array(:); linR_Array_nomut(:)]);
% round up the max value to multiples of 0.5
y_up_lim_LinRich = ceil(maxLinRich/0.5)*0.5;
% y_up_lim_LinRich = 10.5;
maxArrivalR = max([gmeanARank_Array(:); gmeanARank_Array_nomut(:)]);
% round up the max value to multiples of 5
y_up_lim_ArrR = ceil(maxArrivalR/5)*5;
% y_up_lim_ArrR = 45;


figure
subplot_i = 0;
for i_imm_freq=1:length(list_imm_freq)
    subplot_i = subplot_i + 1;
    for i_E=1:length(listE)
        subplot(2,2,subplot_i)
        
        plot(MultiM_opt.Roundedtime,gmeanARank_Array(:,i_E,i_imm_freq),'Color',colorsforplot1{i_E},'LineStyle',Lstyle1{i_E},'LineWidth',2)
        hold on
        plot(MultiM_opt.Roundedtime,gmeanARank_Array_nomut(:,i_E,i_imm_freq),'Color',colorsforplot2{i_E},'LineStyle',Lstyle2{i_E},'LineWidth',1.5)
        
        
        xlim([0 MultiM_opt.Roundedtime(end)]);
        ylim([0 y_up_lim_ArrR]);
        hold on
        ylabel('grand mean arrival rank');
        xlabel('time (in generations)')
        set(gca,'FontSize',14);
        set(gca,'FontName','Ariel')
        set(gca,'box','off','color','none')

        
        
        subplot(2,2,subplot_i+2)

        plot(MultiM_opt.Roundedtime,linR_Array(:,i_E,i_imm_freq),'Color',colorsforplot1{i_E},'LineStyle',Lstyle1{i_E},'LineWidth',2)
        hold on
        plot(MultiM_opt.Roundedtime,linR_Array_nomut(:,i_E,i_imm_freq),'Color',colorsforplot2{i_E},'LineStyle',Lstyle2{i_E},'LineWidth',1.5)
                
        xlim([0 x_up_limit]);
        ylim([0 y_up_lim_LinRich]);
        
        hold on
        ylabel('mean lineage richness')
        xlabel('time (in generations)')
        set(gca,'FontSize',14);
        set(gca,'FontName','Ariel')
        set(gca,'box','off','color','none')
    end
    

end
