function [y, t, Np] = simulation(C0, tspanmax, k)
% Wrapper function for ode solver of the reaction mechanism

% Assign rate constants
k1 = k(1);      k2 = k(2);     k3 = k(3);
k1m = k(4);     k2m = k(5);    k3m = k(6);
k3p = k3 / (5e-7);

% Define timespan
tspan = [0, tspanmax];


% All concentration are nonnegative
opts = odeset('NonNegative',1:10);
% Set tolerance of OSE simulation
opts = odeset(opts,'RelTol',1e-16,'AbsTol',1e-17);
% Run ODE solver
[t, y] = ode15s(@mechanism, tspan, C0, opts);

function dC = mechanism(~, C)
    % Modelling ubiqutine chain reaction

    % User Friendly variable names
    CU = C(1);      CE1 = C(2);     CE1U = C(3);
    CE2 = C(4);     CE2U = C(5);    CE3 = C(6);
    CE3U = C(7);    CE3pU = C(8);   CE1E2 = C(9);
    CE2E3 = C(10);

    % Rate laws
    r1 = k1*CU*CE1;
    r1m = k1m*CE1U.^2;
    r21 = k2*CE1U*CE2;
    r22 = k2*CE1*CE2; 
    r2m = k2m*CE1E2;
    r31 = k3p*CE2U*CE3.^2;
    r32 = k3*CE2U*CE3U;
    r33 = k3*CE2U*CE3pU;
    r34 = k3p*CE2*CE3.^2;
    r3m = k3m*CE2E3;

    % Mass balances
    dCU = -r1 + r1m;
    dCE1 = -r1 + r21 - r22 + r2m + r1m;
    dCE1U = r1 - r21 - r1m;
    dCE2 = -r21 - r22 + r2m + r31 + r32 + r33 - r34 + r3m;
    dCE2U = r21 - r31 - r32 - r33;
    dCE3 = -r31 - r34 + r3m;
    dCE3U = r31 - r32;
    dCE3pU = r32;
    dE1E2 = r22 - r2m;
    dE2E3 = r34 - r3m;

    % Assign output variables
    dC(1,:) = dCU;
    dC(2,:) = dCE1;
    dC(3,:) = dCE1U;
    dC(4,:) = dCE2;
    dC(5,:) = dCE2U;
    dC(6,:) = dCE3;
    dC(7,:) = dCE3U;
    dC(8,:) = dCE3pU;
    dC(9,:) = dE1E2;
    dC(10,:) = dE2E3;

end

% Calculate degree of polymerisation
yE3 = y(:,7) + y(:,8);
yUb = C0(1) - (y(:,1) + y(:,3) + y(:,5));
Np = yUb ./ yE3;

end

