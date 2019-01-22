%% run_rushton_model.m
% Author: Nathan Owen
% Last updated: 14/01/2019

% Clear command line and workspace
clc
clear

% Load clhs experimental design
design = dlmread('clhs_design.txt',' ',1,0);

% Design size is the number of rows in design
design_size = size(design, 1);

% Pre-allocate matrix to store model runs
% Each run returns a time series of length ndays
ndays = 2922;
model_runs = zeros(design_size, ndays);

% Run the Rushton model for each experimental design point in a loop
% Make sure to send through land uses in percentages not proportions
for i = 1:design_size
    model_runs(i,:) = rushton_model(100 * design(i,:));
end

% Save model runs to text file
dlmwrite('model_runs.txt', model_runs, ' ')
