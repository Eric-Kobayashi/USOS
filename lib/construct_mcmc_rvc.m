function data = construct_mcmc_rvc()
% Extract data from data table and construct ydata

% Raw data
addpath('raw_data');
E1_T = readtable('E1_s.csv');
E2_T = readtable('E2_s.csv');
E3_T = readtable('E3_s.csv');

data.C0 = []; data.Edim = [];
for enzyme = [1, 2, 3]
    switch enzyme
        case 1
            dataT = E1_T;
            Vary_E = table2array(dataT(:, 1)); % varying concentration
            Vary_dim = length(Vary_E);
            data.ydata = table2array(dataT(:, 2));
            % Define initial concentrations,   
            C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], Vary_dim, 1);
            % [CU CE1 CE1U CE2 CE2U CE3 CE3U CE3pU CE1E2 CE2E3]
            C0 = C0_mat;
            C0(:, 2) = Vary_E;
        case 2
            dataT = E2_T;
            Vary_E = table2array(dataT(:, 1));
            Vary_dim = length(Vary_E);
            C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], Vary_dim, 1);
            C0 = C0_mat;
            C0(:, 4) = Vary_E;
        case 3
            dataT = E3_T;
            Vary_E = table2array(dataT(:, 1));
            Vary_dim = length(Vary_E);
            C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], Vary_dim, 1);
            C0 = C0_mat;
            C0(:, 6) = Vary_E;
    end
    data.C0 = [data.C0; C0];
    data.Edim = [data.Edim Vary_dim];
    if enzyme ~= 1
        y = table2array(dataT(:, 2));
        y_dim = length(data.ydata);
        a = max(y_dim, Vary_dim);
        data.ydata = [[data.ydata;zeros(a-y_dim, 1)] [y;zeros(a-Vary_dim, 1)]];
    end
end

end

