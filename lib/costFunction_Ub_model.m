function J = costFunction_Ub_model(logk, logtspan_est, E1_T, E2_T, E3_T)
% Calculate the normalised difference between model and experimental data
% For optimising rate constant 

% Transform logk back to k
k = 10.^logk;

% Calculating the costFunctions for each enzyme
J1 = costFunction_enzymes(1, E1_T, k, logtspan_est);
J2 = costFunction_enzymes(2, E2_T, k, logtspan_est);
J3 = costFunction_enzymes(3, E3_T, k, logtspan_est);

J = sum([J1 J2 J3]);

end