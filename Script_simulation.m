close all;
addpath('lib');
% Define simulation time in seconds (Experimental runtime: 6000 sec)
tspanmax = 6000;

% Define rate constants (log10 of the k)
logk0 = [2.9946    7.2742    3.0916   10.7857    2.5059   -0.0490];
k = 10.^(logk0);

% Plot the concentration vs time curves output by the model
% At different enzyme concentration
% Modify the code in plot_model_conc to visualize other species
plot_model_conc(1, k, tspanmax);
% 1 for E1, 2 for E2, 3 for E3



