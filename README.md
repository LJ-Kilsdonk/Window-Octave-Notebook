# Window of opportunity for establishment
This is the code associated with the paper "Transient eco-evolutionary dynamics and the window of opportunity for establishment of immigrants". It was run on Matlab 2015a.  

Two files contain the command-line code to perform the simulations and generate the figures (both main text and appendix):  

- script_for_dual_migration.m- This is the code used for the dual-migration setup  

- script_for_multi_migration.m- This is the code used for the multi-migration setup  

Both files are divided in sections. Each section can be run independently of the rest of the script.  
Because some simulations can take an extensive amount of time we have saved the output of each of these long simulations.  
To redo the simulations in the sections with long simulations, change isredo=false to isredo=true in that section.  


## These are the main functions/scripts to run the model:  

ode_ColonizationDynamics_standard.m- This is a function containing the ODEs of the standard (infinite population size) model. Input is the current state, output is the rate of change for each element of N.   

ode_ColonizationDynamics_stochastic_finite- This is a function containing the stochastic ODEs (finite population size) model. Input is the current state, output is the rate of change for each element of N  

usemodel.m- This is a function used to load the model: for the model provided as input, it makes the parameters, settings, and variables global and gives them their default/initial values (incl. a function handle to the model).  

loadmodel.m- This is command-line code to load the model (specified by the variable modelname) into the workspace. It calls usemodel.m, and makes the parameter/settings/variables global in the current workspace.  

runM.m- This is a function used to run the model (excluding immigration events). It calls the in the settings specified solver, to numerically solve the (stochastic) ODE model over time.  

runMgrow.m- This calls the runM.m and adds one column to N, to allow for the addition of an additional immigrant.  

dualmigrunandplot.m- This is a function that does a run of the model with two immigration events (dual-migration setup) and produces a muller plot. It uses runM.m to run the model inbetween and after the two immigration events.  

multi_migration_run.m-  This is a function that does a run of the model with many immigration events (multiple-migration setup). It calls runMgrow to add a new column to N before each immigration event. By doing so, for most of the run, N has few columns, which considerably speeds up the solver. If there are no mutations, each lineage has only one phenotype. To speed up these simulations, this function turns N into a column vector of lineages, each with 1 phenotype and death rate.  



## These are the functions that were used to explore parameter space:  

WindowSearch.m- Finds the window of opportunity for establishment for each combination of founder and invader phenotype using the build-in function fzero.m. It returns a matrix with the time window for each combination of phenotypes.  
    
* plot_window.m- Plots the window of opportunity using the matrix returned by WindowSearch.m  

adaptation_popsize.m- Tests the time needed to adapt with the finite population size model for different population sizes  

loopmulti.m- Loops over different settings (E-values and immigration frequencies) and replicates, each iteration it calls multi_migration_run and stores the final population density of each lineage.  

 * linrichness.m- Creates an errorbar plot of the lineage richness using the array of final population densities per lineage at different settings as returned by loopmulti.m  

* arrivalrank.m- Creates a stairplot of the mean arrival rank using the array of final population densities per lineage at different settings as returned by loopmulti.m. It also displays the grand mean arrival rank for runs with mutations and the difference in grand mean arrival rank with compared to without mutations (i.e. the effect size).  

preequilloopmulti.m- A loop like loopmulti.m, but instead of storing only the final values, it stores the change over time. To reduce memory usage, after each set of replicates (with the same settings) it immediately calculates a grand mean arrival rank and mean lineage richness. Unlike loopmulti.m, which has to be called twice: once with and once without mutations, this function does run with and without simulations each loop.  
    
* preEquilFigs.m- Creates plots over time of the mean lineage richness and grand mean arrival rank, using the output of preequilloopmulti.m.  



## These are the minor functions used to run the model:    

updateQandd.m- This function updates the mutation matrix Q, and the death rates d. For the multiple-migration setup, if there are no mutations it adapts Q and d for the column vector version of N. It is affected by changes in many of the settings, and thus best called before each run.  

enterfounder.m- This function introduces the first immigrant by changing state variable N to a matrix of zeros and adding one propagule density to the specified phenotype of the first lineage.  

change_dimsN.m- This function updates the dimensions of N (and the associated settings in the model).  

store.m- This is a function that can be used to keep the final value after a run. runM.m first stores the output of the run over time in my_N. This function changes the state variable N to the final value at the end of a run.  

euler.m- Solves ODEs using the Euler method with fixed time steps. Used for numerically solving the stochastic finite version over time.  

n1maxtest.m- This function tests for the maximum number of lineages in a simulation. For memory allocation in loopmulti.m and preequilloopmulti.m, a maximum number of lineage has to be set. Randomly drawing arrival times over, a distribution of maximum arrival ranks is calculated and plotted.  

getknownfield.m- This function tests if field f of a structure exists, and if not, gives it value provided as input argument. This function is used to assign default values to the fields of structures.  


