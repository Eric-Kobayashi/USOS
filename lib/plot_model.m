function [Test_E, y, Vary_E, y0] = plot_model(enzyme, k, logtspan_est)
% Plot model vs experimental Ubrate
% Import raw data for plotting
addpath('raw_data');
E1_T = readtable('E1.csv');
E2_T = readtable('E2.csv');
E3_T = readtable('E3.csv');

Test_E = linspace(-7.5, -4, 100);
Test_E = 10.^Test_E;

% Changing E
switch enzyme
    case 1
        dataT = E1_T;
        Vary_E = table2array(dataT(:, 1));
        % Define initial concentrations   
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Test_E), 1);
        % [CU CE1 CE1U CE2 CE2U CE3 CE3U CE3pU CE1E2 CE2E3]
        C0 = C0_mat;
        C0(:, 2) = Test_E;
    case 2
        dataT = E2_T;
        Vary_E = table2array(dataT(:, 1));
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Test_E), 1);
        C0 = C0_mat;
        C0(:, 4) = Test_E;
    case 3
        dataT = E3_T;
        Vary_E = table2array(dataT(:, 1));
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Test_E), 1);
        C0 = C0_mat;
        C0(:, 6) = Test_E;
end

y0 = table2array(dataT(:, 2)); % measured Ub rate
Vary_dim = size(C0, 1);
y = zeros(Vary_dim, 1);
tspan = 10.^logtspan_est;

for i = 1:Vary_dim
    [Ubc, tpoint, ~] = simulation(C0(i, :), tspan, k);
    [y(i), ~, ~] = calcUbrate(Ubc, tpoint);
end

end

