function plot_mcmc(enzyme, k_arr, logtspan_est)

addpath('raw_data');
E1_T = readtable('E1.csv');
E2_T = readtable('E2.csv');
E3_T = readtable('E3.csv');
Ub_T = readtable('Ub.csv');

Test_E = linspace(-7.5, -4, 100);
Test_E = 10.^Test_E;

% Changing E
switch enzyme
    case 1
        dataT = E1_T;
        Vary_E = table2array(dataT(:, 1));
        % Define initial concentrations,   
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Test_E), 1);
        % [CU CE1 CE1U CE2 CE2U CE3 CE3U CE3pU]
        C0 = C0_mat;
        C0(:, 2) = Test_E;
        enzyme_label = 'E1';
    case 2
        dataT = E2_T;
        Vary_E = table2array(dataT(:, 1));
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Test_E), 1);
        C0 = C0_mat;
        C0(:, 4) = Test_E;
        enzyme_label = 'E2';
    case 3
        dataT = E3_T;
        Vary_E = table2array(dataT(:, 1));
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Test_E), 1);
        C0 = C0_mat;
        C0(:, 6) = Test_E;
        enzyme_label = 'E3';
    case 0
        dataT = Ub_T;
        Vary_E = table2array(dataT(:, 1));
        Test_E = linspace(-6, -4, 100);
        Test_E = 10.^Test_E;
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Test_E), 1);
        C0 = C0_mat;
        C0(:, 1) = Test_E;
        enzyme_label = 'Ub';
end

y0 = table2array(dataT(:, 2)); % measured Ub rate
Vary_dim = size(C0, 1);
y = zeros(Vary_dim, 1);
tspan = 10.^logtspan_est;

disp(enzyme_label);
figure;  
% C = 1;
for n = 1:size(k_arr, 1)
    k = k_arr(n, :);
    for i = 1:Vary_dim
        [Ubc, tpoint, ~] = simulation(C0(i, :), tspan, k);
        [y(i), ~, ~] = calcUbrate(Ubc, tpoint);
    end

    % Normalise y to y0
    ave_y = mean(y);
    ave_y0 = mean(y0);
    C = ave_y0/ave_y;
    y_norm = C * y;

    slx = semilogx(Test_E, y_norm(:,1), ':');
    alpha(slx, .2);
    hold on;
end

scatter(Vary_E, y0(:,1));
title(strcat('Vary ', enzyme_label, ' concentration')); hold off;

% figure;
% semilogx(Test_E, tchoice);


end
