clc; clear; close all;
fp1 = fopen('.\data.txt', 'w');
%% Generate connection matrix============================
TT = 80;
K = 2;
for i = 1:21
    Per = 0.05 * (i-1);
    %% Run multiple simulations (20000 different random seeds)
    num_runs = 20000;
    seeds = randi([1, 100000], num_runs, 1);  % Generate 20000 random integer seeds between 1 and 100000
    degree_values = zeros(num_runs, 16);

    for run = 1:num_runs
        rng(seeds(run));  % Set random seed for current run
        [matrix, x, y] = func_WS_network(TT, K, Per);
        [Dds, Dds_avg, M, P_Dds] = func_Degree_Distribution(matrix);
        degree_values(run, 1:length(P_Dds)) = P_Dds;
    end

    degree_average = sum(degree_values, 1) / num_runs;

    % Input data
    degree = 0:15;  % Numerical column

    % Calculate mean μ
    mu = sum(degree .* degree_average);

    % Calculate variance σ²
    sigma_squared = sum((degree - mu).^2 .* degree_average);

    fprintf(fp1, '%f %f\n', Per, sigma_squared);
end

fclose(fp1);
