% Set different prior for the rate constant
warning off;
rng('shuffle');
addpath('lib');
% Define initial k
Orig_k = [2.9946    7.2742    3.0916   10.7857    2.5059   -0.0490];
%         k1        k2        k3       k1m        k2m      k3m
dim_k = length(Orig_k);

% Number of epochs
epoch = 4;

% Define simulation time in seconds
logtspan_est = log10(6000);

% Read experimental data
addpath('raw_data');
E1_T = readtable('E1.csv');
E2_T = readtable('E2.csv');
E3_T = readtable('E3.csv');

logk_min_arr = zeros(epoch, 4);
J_arr = zeros(epoch, 1);
% Each iteration has a slightly different initial state
deltak = 0.5*rand(epoch-1, dim_k) - 0.25;
Prior_k = [Orig_k; Orig_k+deltak];

% Use multi-core process for simulated annealing
parfor i = 1:epoch
    J_init = costFunction_Ub_model(Prior_k(i, :), logtspan_est, E1_T, E2_T, E3_T)
    try
        [logk_min_arr(i, :), J_arr(i), output] = Vary_all_E_model(Prior_k(i, :), logtspan_est);
    catch
        warning('Error in optimization. Skip the loop');
        logk_min_arr(i, :) = Prior_k(i, :);
        J_arr(i) = J_init;
        continue;
    end
    disp(["Iters: ", output.iterations, "J = ", J_arr(i)]);
end
