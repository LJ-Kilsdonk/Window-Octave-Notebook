function dualmigrunandplot(colonizer_pheno,invader_pheno,lag_invasion,time_lengthrun,varargin)
%DUALMIGRUNANDPLOT  runs with two migration events with specified phenotype
%   and time lag between events (in generations).
%   Requires the model to be loaded with function usemodel or command-line
%   script loadmodel.
%   
%   DUALMIGRUNANDPLOT(colonizer_pheno,invader_pheno,lag_invasion,time_lengthrun,)
%   where colonizer_pheno is the phenotype of the first migrant,
%   invader_pheno the phenotype of the second migration, lag_invasion the
%   lag in time between the arrival of the two migrant, time_lengthrun the
%   duration of (in generations) of the run
% 
%   DUALMIGRUNANDPLOT(...,'PropertyName',PropertyValue,...) sets the value
%   of the specified plot property. The two properties are colMAP (by
%   default jet), which is the colormap for the muller plot and
%   whiteSpaceColor (by default kind of white: [0.886,0.902,0.914]) which
%   sets the  background color of the muller plot. 
%   
%   See also: loadmodel, usemodel, runM, multi_migration_run
%
global propsize my_N my_t N a


if ~exist('propsize','var')||isempty(propsize)
    propsize=0.001;
    disp('no propagule_size defined, propagule_size set to 0.001')
end     
% time steps plotted (note internally time steps can differ).
tstep_size = 0.1;

% Set all elements N to zero and introduce first migrant.
enterfounder(colonizer_pheno,propsize) 
% Run model till the arrival of the invader or time_lengthrun (whichever
% comes 1st).
runM(0:tstep_size:min(time_lengthrun,lag_invasion),1) 
t_preInv = my_t; % time data of this period before the arrival of the invader
N_preInv = my_N; % data N (over time) of this period before the arrival of the invader

% Introduce invader
N(invader_pheno,2) = propsize;

% Run model till end of simulation:
if (time_lengthrun - lag_invasion) > 0
    if (time_lengthrun - lag_invasion) > tstep_size
        runM(lag_invasion:tstep_size:time_lengthrun)
    elseif (time_lengthrun - lag_invasion) > 0
        runM([lag_invasion time_lengthrun])
    end
    t_postInv = my_t(2:end); % time vector of this period after (EXCLUDING THE MOMENT ITSELF) the arrival of the invader
    N_postInv = my_N(2:end,:); % data N (over time) of this period after the arrival of the invader
    % Combine both runs
    my_t = [t_preInv;t_postInv];
    my_N = [N_preInv;N_postInv];
end

% Default settings for the mullerplot
if sum(strcmp(varargin,'colMAP'))==1
    i_var=find(strcmp(varargin,'colMAP'));
    MAP = varargin{i_var+1};
else
    MAP = jet(a);
end
if sum(strcmp(varargin,'whiteSpaceColor'))==1
    i_var=find(strcmp(varargin,'whiteSpaceColor'));
    whiteSpaceColor = varargin{i_var+1};
else
    whiteSpaceColor = [0.886,0.902,0.914];
end

% plot results run as muller plot:
mullerplot([1 0 1],MAP,whiteSpaceColor)


%------------------Muller plot--------------------------

function mullerplot(Wratio,MAP,whiteSpaceColor,varargin)
%MULLERPLOT    Produce muller plots of a run with the dual migration setup
% 
%   MULLERPLOT(Wratio, MAP, whiteSpaceColor, 'Name',Property)
%   where Wratio is a vector of size 3 with the ratio of the blank space
%   above, in between, and below both lineages. MAP is the colormap used,
%   whiteSpaceColor the color of the blank space, and a name property pair
%   can be used to change additional settings from deafult (i.e. the
%   resolution). The default resolution is 1024.

global a my_N my_t
% Default settings
if sum(strcmp(varargin,'resolution'))==1
    indvar=find(strcmp(varargin,'resolution'));
    resolution=varargin{indvar+1};
else
    resolution=1024;
end
if length(Wratio)~=3
    error('Wrat needs 0 or 3 input arguments, but was given %d',length(varargin{1}))
else
    var1=Wratio;
    Wrat=var1/sum(var1); % ratio of white space (adding up to 1).
end

% my_N is a matrix, with rows being different points in time during the run
% and the columns describing both the different phenotypes and lineages
% (they are the elements of N ordered like N(:).
leafmut = my_N(:,1:a); % densities by phenotype for lineage 1 (time x phenotype)
leafmut2 = my_N(:,a+1:end); % densities by phenotype for lineage 2 (time x phenotype)
maxWidth = max(sum(my_N,2)); % the maximum total population density, is the width (difference between the min and max ylim).
WSPACE=(maxWidth-sum(my_N,2)); % blank space over time
numberofcolumns = 4*a +3; % 2 lineages x 2 parts (aka leaf) in plot (top and bottom part/leaf) plus 1 blank space below, 1 blank space above and one blank space in between both lineages
edgeColumns = (numberofcolumns-1)/2;
edgeRGB = [1 1 1]; 

% Plot
barsmooth(my_t,...
    [WSPACE*Wrat(1) leafmut(:,end:-1:1)/2 leafmut/2 WSPACE*Wrat(2) leafmut2(:,end:-1:1)/2 leafmut2/2 WSPACE*Wrat(3)],...
    [whiteSpaceColor; MAP(end:-1:1,:); MAP; [1 1 1]; MAP(end:-1:1,:); MAP;whiteSpaceColor],...
    resolution,  edgeColumns,edgeRGB)

% Adjust plot limits and add labels
ylim([0 maxWidth]);
xlim([0 my_t(end)])
xlabel('time')
ylabel('population density')
set(gca,'FontSize',14);
set(gca,'FontName','Ariel')

% add seperation line between the two lineages
start1stline = find(sum(leafmut2~=0,2),1);
hold on
% Plot first a white line, and then a dashed line on top. White line makes
% the dashed line better visible
plot(my_t(start1stline:end),(WSPACE(start1stline:end)*Wrat(1)) + sum(leafmut(start1stline:end,:),2),'-w','LineWidth',1)
plot(my_t(start1stline:end),(WSPACE(start1stline:end)*Wrat(1)) + sum(leafmut(start1stline:end,:),2),'--k','LineWidth',1)
set(gca,'YtickLabel',[]); % Remove vertical axis labels


%-----------------------------barsmooth--------------------------------
function barsmooth(X,Y,CMAP,resolution,edgeColumns,edgeRGB)
%BARSMOOTH  Creates a smooth stacked bar plot with no gaps. Smooth in the
%   sense that, when the locations of the different bars are not equidistant, this
%   function unlike the bar plot, does make them fit exactly together. (each
%   bar is made out the region between two points, rather than being the
%   middle of one. Also the bar is not flat. it has an inclined top. Such that
%   actually the sides to correspond to the actual values at those points in
%   time and in between there is a linear approximation.

% setting the maximum resolution (so that very fine detailed parts get
% erased, because you can't see it on the screen anyway)
if nargin<4
    resolution=1024;
end
if nargin<6
    edgeColumns=[];
    edgeRGB=[];
end

% Reducing the resolution if it is too high
toofine=0;
i=0;
while ~isempty(toofine)
    toofine = find(X(3:2:end)-X(1:2:end-2)<((max(X)-min(X))/resolution))*2;
    X(toofine)=[];
    Y(toofine,:)=[];
    i=i+1;
    if i>1000
        warning('resolution was not reduced properly due to iterations reaching a limit')
        break
    end
end

Yc_h=cumsum(Y,2);
Yc_l=[zeros(size(Y,1),1) Yc_h(:,1:end-1)];

hold on

Faces=1:length(X)*2;
for j=1:size(Yc_h,2)
    Verts=[[X; X(end:-1:1)], [Yc_l(:,j); Yc_h(end:-1:1,j)]];
    if mean( Yc_h(1:end,j) - Yc_l(:,j) <= .1/resolution)~=1       
        if sum(j==edgeColumns)
            patch('Faces',Faces,'Vertices',Verts,'EdgeColor',edgeRGB(j==edgeColumns,:),'FaceColor','flat','FaceVertexCData',CMAP(j,:));
        else
            patch('Faces',Faces,'Vertices',Verts,'EdgeColor','none','FaceColor','flat','FaceVertexCData',CMAP(j,:));
        end
    end
end
hold off
