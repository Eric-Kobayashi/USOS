function [logk_min, J, output] = Vary_all_E_model( logk0, logtspan_est)

% Read experimental data
addpath('raw_data');
E1_T = readtable('E1.csv');
E2_T = readtable('E2.csv');
E3_T = readtable('E3.csv');

% set cost function
Jfunc = @(logk)(costFunction_Ub_model(logk, logtspan_est, E1_T, E2_T, E3_T));

% Use simmulated annealing for optimizing rate constants
% Lower and higher bound for estimating rate constants
lb = [2.5   0  0   7  -4  -4];
ub = [4.5  10  8  14   8   8];
% Initial temperature for simulated annealing
options = optimoptions(@simulannealbnd,'InitialTemperature',25);
[logk_min, J, ~, output] = simulannealbnd(Jfunc, logk0, lb, ub, options);

end

