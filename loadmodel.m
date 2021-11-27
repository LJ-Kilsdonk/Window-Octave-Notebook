%% Load model
% makes the parameters, variables, function handles and settings of the
% model global, and allow access to these variables from the command line.
%
% Before calling loadmodel, the variable modelname has to be defined!
%
global r mut n1 a E d d_min d_diff Q K h propsize   % parameters
global MY_SETTINGS                                  % settings
global N                                            % state variables
% Load packages for Octave
pkg load statistics

if exist('modelname','var')
    usemodel(modelname);
    disp(['loaded model: ',MY_SETTINGS.modelname]);
else
    error('The variable modelname was not defined before calling loadmodel.m')
end