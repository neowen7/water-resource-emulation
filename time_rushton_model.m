%% time_rushton_model.m
% Author: Nathan Owen
% Last updated: 14/01/2019

% Clear command line and workspace
clc
clear

% Load clhs experimental design
design = dlmread('clhs_design.txt',' ',1,0);

% Set position parameter 'pos' = 6, to indicate Frome basin
pos = 6;

% Design size is the number of rows in design
design_size = size(design, 1);

% Set number of timing repetitions
no_of_timings = 100;

% Pre-allocate matrix to store model run times for each repetition
timing_results = zeros(no_of_timings,1);

% Run the experiment no_of_timing times, record run time for each
for timing = 1:no_of_timings
    tic
    for i = 1:design_size
        model05092017fun(100 * design(i,:), pos);
    end
    timing_results(timing) = toc;
end

% Take average of model run times
mean(timing_results)

% Model average run time: 8.7684
% To be compared with emulator average run time

