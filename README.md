---
title: USOS - Ubiquitin System ODE Simulator
author: Eric Hidari
date: 30 December 2019
---

# USOS
MATLAB model for ubiquitin chain reaction detailed in the paper:
A general _in vitro_ assay to study enzymatic activities of the ubiquitin system in __Biochemistry__ USOS uses machine learning method to infer the rate constants in the ODE functions given the FRET measurement data. Namely, simulated annealing, an advanced Markov Chain Monte Carlo (MCMC) algorithm is used here. 

# Prerequisite
- MATLAB R2018b or higher
- MATLAB Global Optimization Toolbox, Parallel Computing Toolbox and Curve fitting Toolbox

# Usage
## Plot initial reaction rate vs E1/E2/E3 concentration
Run Script_plot_model.m with the defined rate constants.

## Plot species concentration vs time at different E1/E2/E3 concentration
Run Script_simulation.m with the defined rate constants. Modify the code in lib/plot_model_conc.m to simulate concnetrations of different species.

## Run simulated annealing to optimize rate constants
Run Script_optimise_rate_constants.m to infer rate constants from the data defined in raw_data folder. Change the epoch number in the script to run from multiple initial states.
The raw_data csv file is arranged as follows:
- The first column is the enzyme concentrations, 
- The second column is the measured initial reaction rates, 
- The third column is the standard deviation of the measured rates (unused).
The optimization algorithm simulated annealing can accept lower/upper bounds of the rate constants and the initial temperature as the parameters. These can be modified in Vary_all_E_model.m

## Run MCMC to estimate rate constants posterior distribution
Alternatively, one can run Metropolis MCMC algorithm to observe the posterior distribution of the rate constants given the data in a Bayesian approach. Run Script_mcmc.m and change the chain length (nsimu) to generate a large number of posterior observations. 

# Dependency
The MCMC package is written by Marko Laine and is downloaded from:
https://mjlaine.github.io/mcmcstat/

# License
This software package is under the MIT License. See LICENSE.txt for details.

