function plot_window(windowOfOpportunity,colorforlines,lineIndeces2plot)
%PLOT_WINDOW  plots window of opportunity for establishment
%   Requires a matrix with the window of opportunity for establishment
%   windowOfOpportunity, with rows corresponding to different phenotypes
%   for the founder and columns for different phenotypes of the invader
%   (i.e. 2nd migrant). Returns a plot of the window of opportunity for
%   establshment. lineIndeces2plot is used to only plot the window of
%   opportunity for invasion with a selection of founder phenotypes,
%   instead of all possible ones.
%   
%   PLOT_WINDOW(windowOfOpportunity,colorforlines,lineIndeces2plot) plots
%   the windowOfOpportunity (the output from WindowSearch with the founder
%   phenotypes corresponding to lineIndeces2plot and using colorforlines
%   (by default jet(a)) as color to plot the window of opportunity (as a
%   line) for each of these founder phenotypes. By default all lines are
%   plotted. 
%
%   See also: WindowSearch
%

%
global d a d_min
% In plot, each column is a line. To plot a window of opportunity line for
% each phenotype of the founder, windowOfOpportunity' is transposed.
window = windowOfOpportunity';
% default settings:
if nargin < 3
    n_lines = a;
    lineIndeces2plot = 1: n_lines;
    if nargin < 2
        colorforlines = jet(a);
    end
else
    n_lines = length(lineIndeces2plot);
end


if d_min == 0
    y_axVals = repmat(-d(:,1),1,n_lines);
%     my_ylabel = 'phenotype-environment mismatch of invading immigrant';
%     my_ylabel_colorbar = 'phenotype-environment mismatch of founder';
    my_ylabel = 'resource independent death rate of invading immigrant';
    my_ylabel_colorbar = 'resource-independent death rate of founder';
    
else
    y_axVals = repmat(-d(:,1) + d_min,1,n_lines);
    my_ylabel = 'resource independent death rate of invading immigrant (in difference from optimum)';
    my_ylabel_colorbar = 'resource-independent death rate of founder (in difference from optimum)';
end


% window contains '0' wherever there is no window of opportunity for
% establishment. Given the horizontal axis is the window of opportunity,
% two or more consecutive 0 values create a vertical line on top of the
% vertical axis. The figure is clearer without these vertical lines.
% Therefore, we turn all the 'zeros' to  NaN, except the first zero (i.e.
% the 'no' establishment value belonging to the most pre-adapted phenotype
% that is never able to invade). 
for iPhF = 1 : a    % iPhF index phenotype founder (1st immigrant)
    i_firstZero = find(window(:,iPhF)==0,1,'first');
    window(window(:,iPhF)==0,iPhF) = NaN; % all 0 values to NaN
    window(i_firstZero,iPhF) = 0; % restore first zero value to 0.
end


f_handle =figure;
plot1 = plot(window(:,lineIndeces2plot), y_axVals);
% Color the lines:
for i=1:n_lines
    set(plot1(i),'Color',colorforlines(i,:))
end

title('The window of opportunity for establishment','FontSize',14)

my_colorbar = colorbar('Location','East');
colormap(colorforlines(end:-1:1,:))
xlabel('generations','FontSize',14)
ylabel(my_ylabel,'FontSize',14)
ylabel(my_colorbar,my_ylabel_colorbar,'FontSize',14)
ylim([-1 0])
% set(f_handle.CurrentAxes,'YTickLabel',sprintfc('%0.1f',-get(f_handle.CurrentAxes,'YTick')'))
% set(my_colorbar,'TickLabels',my_colorbar.TickLabels(end:-1:1))
