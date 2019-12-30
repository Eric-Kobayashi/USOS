close all;
warning off;
addpath('lib');
% Degree of Ub polymerisation

% Define rate constants (log10 of the k)
logk_min = [2.9946    7.2742    3.0916   10.7857    2.5059   -0.0490];
%           k1        k2        k3       k1m        k2m      k3m
% k3prime is automatically derived from k3
k = 10.^ logk_min;

% Define simulation time in seconds (Experimental runtime: 6000 sec)
logtspan_est = log10(6000);

% Raw data from FRET measurement
addpath('raw_data');
E1_T = readtable('E1.csv');
E2_T = readtable('E2.csv');
E3_T = readtable('E3.csv');

% Calculate costFunction
J = costFunction_Ub_model(logk_min, logtspan_est, E1_T, E2_T, E3_T);
fprintf("Mean Square Error of fitting: %s", num2str(J));

e = [3 2 1];

% x, y - simulation; x0, y0 - data
x = {}; y = {}; x0 = {}; y0 = {};

% Simulate the reactions
parfor i=1:length(e)
    warning off;
    [x{i}, y{i}, x0{i}, y0{i}] = plot_model(e(i), k, logtspan_est);
end

% Plot model vs data
figure;
for i=1:length(e)
    subplot(1,3,i);
    x_plot = x{i}*1e6;
    y_plot = smooth(y{i}*1e6);
    semilogx(x_plot, y_plot, ':', 'LineWidth',0.8,'color','black'); hold on;

    switch e(i)
        case 1
            enzyme_label = 'Uba1';
            err = table2array(E1_T(:, 3));
            cl = 'blue';
        case 2
            enzyme_label = 'UbcH5';
            err = table2array(E2_T(:, 3));
            cl = [0 0.7 0.7];
        case 3
            enzyme_label = 'CHIP';
            err = table2array(E3_T(:, 3));
            cl = [0 0.7 0];
        case 0
            enzyme_label = 'Ub';
            err = table2array(Ub_T(:, 3));
            cl = 'blue';
    end 
    errorbar(x0{i}*1e6, y0{i}*1e6, err*1e6, 's',...
    'color', cl, 'MarkerFaceColor',cl,...
    'MarkerSize',4);
    ylim([0 0.0020]);
    xticks([0.01 0.1 1 10 100]);
    xticklabels([0.01 0.1 1 10 100]);
    ylabel('Initial rate / µM s^{-1}');
    xlabel(strcat('[', enzyme_label, '] / µM'));
    title({'Ubiquitination rates with increasing', strcat('[', enzyme_label, '] and model fitting')}); hold off;
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.8;
    set(ax,'box','off');
    set(gca,'TickDir','out');
end

set(gcf,'position',[161,56,1265,329]);