function final_Np = plot_model_conc(enzyme, k, tspan)
% Plot model vs experimental Ub_conc at a given tspan

% Import raw data for plotting
addpath('raw_data');
E1_T = readtable('E1.csv');
E2_T = readtable('E2.csv');
E3_T = readtable('E3.csv');

% Changing E
switch enzyme
    case 1
        dataT = E1_T;
        Vary_E = table2array(dataT(:, 1));
        % Define initial concentrations
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Vary_E), 1);
        % [CU CE1 CE1U CE2 CE2U CE3 CE3U CE3pU CE1E2 CE2E3]
        C0 = C0_mat;
        C0(:, 2) = Vary_E;
        enzyme_label = 'E1';
    case 2
        dataT = E2_T;
        Vary_E = table2array(dataT(:, 1));
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Vary_E), 1);
        C0 = C0_mat;
        C0(:, 4) = Vary_E;
        enzyme_label = 'E2';
    case 3
        dataT = E3_T;
        Vary_E = table2array(dataT(:, 1));
        C0_mat = repmat([5e-6 5e-7 0 5e-7 0 5e-7 0 0 0 0], length(Vary_E), 1);
        C0 = C0_mat;
        C0(:, 6) = Vary_E;
        enzyme_label = 'E3';
end

Vary_dim = size(C0, 1);
final_Np = zeros(Vary_dim, 1);
cmap = jet(Vary_dim);
for i = 1:Vary_dim
    hold on;
    [y, t, Np] = simulation(C0(i, :), tspan, k);
    figure(1);
    % Change y(:,1) to other species for plotting
    %  1  2   3    4   5    6   7    8     9     10   
    % [CU CE1 CE1U CE2 CE2U CE3 CE3U CE3pU CE1E2 CE2E3]
    plot(t/60, y(:,1)*1e6, ':', 'LineWidth',1.8, 'color', cmap(i, :));
    final_Np(i) = Np(end);
    hold on;
    figure(2);
    plot(t/60, Np, ':', 'LineWidth',1.8, 'color', cmap(i, :));
    % Make legends
    if i == 1
        str = {strcat(enzyme_label, ': ', num2str(Vary_E(i)), ' M')};
        str = [str , strcat(enzyme_label, ': ', num2str(Vary_E(i)), ' M')];
    else
        str = [str , strcat(enzyme_label, ': ', num2str(Vary_E(i)), ' M')];
        str = [str , strcat(enzyme_label, ': ', num2str(Vary_E(i)), ' M')];
    end

end

hold off;

% plot your data
figure(1);
legend(str{:});
xlim([0 tspan/60]);
xlabel('Time (min)');
ylabel('Concentration / µM');
title({'Species concentration vs time with increasing', strcat('[', enzyme_label, ']', ' predicted by model')}); 
set(gcf,'position',[161,596,565,399]);

figure(2);
legend(str{1:2:end});
xlim([0.1 tspan/60]);
xlabel('Time (min)');
ylabel('Degree of polymerization')
title({'Ubiquitination vs time with increasing', strcat('[', enzyme_label, ']', ' predicted by model')}); 
set(gcf,'position',[161,56,565,399]);

end

