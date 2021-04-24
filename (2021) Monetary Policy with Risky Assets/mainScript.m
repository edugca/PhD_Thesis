%%%%%%%%% Amaral(2021): Monetary Policy with Risky Assets

addpath(genpath('Auxiliary functions'));

addpath('Data');
addpath('Filters');
addpath('Power');
addpath('Model');
addpath('Simulation');
addpath('Empirical Correlations');
addpath('Taylor Coefficient');

%% Plot figures 1, 21, and 22
rN_filters;

%% Figure 2

% It is a scheme drawn with a Mac software called Graph Sketcher
% https://github.com/graphsketcher/GraphSketcher
%
% In the folder "Active vs Passive"
% Editable file: Monetary and Fiscal Policy - Active vs Passive.ograph
% Published file: Monetary and Fiscal Policy - Active vs Passive.png

%% Plot figures 3, 4, 25, and 26
partialEquilibrium;

%% Plot figures 23, 24, 27, and 28
partialEquilibrium_riskyNatural;

%% Plot figures 5, 6, 29, and 30
partialEquilibrium_riskyPolicy;

%% Plot figures 7, 8, 9, and 10
powerOfMonPol;

%% Plot figures 11, 12, 36 to 43
simulation;

%% Plot figures 13 to 19
empiricalCorrelations;

%% Plot figure 20
taylorCoeff_infVol;

%% Plot figures 31 to 35
powerOfMonPol_coeff;
