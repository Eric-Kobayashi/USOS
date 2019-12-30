function J = costFunction_enzymes(enzyme, dataT, k, logtspan_est)
%costFunction_enzymes calculates the error between model and experiment
% for one enzyme varying concentration model
Vary_E = table2array(dataT(:, 1)); % varying concentration
y0 = table2array(dataT(:, 2)); % measured Ub rate
Vary_dim = length(Vary_E);

% Define initial concentrations,   
C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], Vary_dim, 1);
% [CU CE1 CE1U CE2 CE2U CE3 CE3U CE3pU CE2E3]
C0 = C0_mat;

% Changing E
switch enzyme
    case 1
        C0(:, 2) = Vary_E;
    case 2
        C0(:, 4) = Vary_E;
    case 3
        C0(:, 6) = Vary_E;
end

tspan = 10.^logtspan_est;

y = zeros(Vary_dim, 1);
for i = 1:Vary_dim
    [Ubc, tpoint, ~] = simulation(C0(i, :), tspan, k);
    [y(i), ~, ~] = calcUbrate(Ubc, tpoint);
end

% Mean square error
J = 1/Vary_dim*sum((y - y0).^2)*1E19;

end

