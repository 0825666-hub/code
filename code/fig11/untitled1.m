clc; clear; close all;
fp1 = fopen('.\data.txt', 'w');
%% Generate connection matrix============================
TT = 80;
K = 2;
Per = 0.85;
%% Run multiple simulations (10000 different random seeds)
num_runs = 10000;
seeds = randi([1, 100000], num_runs, 1);  % Generate 10000 random integer seeds between 1 and 100000
degree_values = zeros(num_runs, 15);

for run = 1:num_runs
    rng(seeds(run));  % Set random seed for current run
    [matrix, x, y] = func_WS_network(TT, K, Per);
    [Dds, Dds_avg, M, P_Dds] = func_Degree_Distribution(matrix);
    degree_values(run, 1:length(P_Dds)) = P_Dds;
end

degree_average = sum(degree_values, 1) / num_runs;

for i = 1:length(degree_average)
    fprintf(fp1, '%f %f\n', i-1, degree_average(i));
end

fclose(fp1);

%% Plot================================
% Plot bar chart
bar(1:15, degree_average, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k', 'BarWidth', 0.8);

% Set axis labels
xlabel('Degree (k)');
ylabel('Proportion of Nodes');
title('Degree Distribution (Pre-computed Proportions)');

% Fix x-axis ticks at each bar center
xticks(1:15);  
grid on;