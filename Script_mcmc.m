% Set up Markov chain Monto Carlo method for Bayesian inference
warning off;
addpath('mcmcstat');
addpath('lib');
clear data model params options

% Prior rate constants
logk0 = [2.9946    7.2742    3.0916   10.7857    2.5059   -0.0490];

% Read concentration vs time data from the experimental measurement
data = construct_mcmc_rvc();

ss0 = costFunction_rvc(logk0, data);
mse = ss0/(length(data.ydata));
fprintf("Mean Square Error of prior parameters: %s", num2str(mse));

% Run parameters
method = 'dram';
nsimu = 2000;
adaptint = 200;

params = {
%   name,  init,        min, max, mu,       sig, target?
   {'k1',  logk0(1),    2.5, 4.5, logk0(1), 1.0, 1       }
   {'k2',  logk0(2),    0,   10,  logk0(2), 2.0, 1       }
   {'k3',  logk0(3),    0,   8,   logk0(3), 2.0, 1       }
   {'k1m', logk0(4),    7,   14,  logk0(4), 2.0, 1       }
   {'k2m', logk0(5),    -4,  8,   logk0(5), 2.0, 1       }
   {'k3m', logk0(6),    -4,  8,   logk0(6), 2.0, 1       }
};

model.ssfun = @costFunction_rvc;
model.sigma2 = mse;

options.method      = method;        % adaptation method (mh,am,dr,dram)
options.nsimu       = nsimu;         % n:o of simulations
options.adaptint    = adaptint; % adaptation interval
options.printint    = 200; % how often to show info on acceptance ratios
options.verbosity   = 1;  % how much to show output in Matlab window
options.waitbar     = 1;  % show garphical waitbar
options.updatesigma = 1;  % update error variance
options.stats       = 1;  % save extra statistics in results

[results,chain,s2chain] = mcmcrun(model,data,params,options);
%save()

