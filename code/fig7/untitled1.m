clc; clear; close all;

%% Parameter Settings
TT = 1;
cm = 1.0;
pi = 3.1415926;

step = 0.01;
Tao = 2;

A = 0.5;          % A cannot be too large
T = 23;           % Temperature factor temporarily not added

m_period = 600;  % Select 600 periods

gc = 0.2;         % Conductance between two compartments
p_ = 0.5;         % Ratio of dendritic area to total cell area

fp1 = fopen('.\data.txt', 'w');

% w value range and settings
w_values = 0.01:0.0022:0.12;  % 50 data points
num_w = length(w_values);
% D value range and settings
D_values = 0:0.8:40;  % 50 data points
num_D = length(D_values);

% Run multiple simulations (10 different random seeds)
num_runs = 10;  % Run 10 times
seeds = [1,2,3,4,5,6,7,8,9,10];  % 10 different random seeds

%% Main loop: iterate through all w values
for w_idx = 1:num_w
    w = w_values(w_idx);
    fprintf('\n=== Current w value: %.4f ===\n', w);
    
    % Calculate T1 and T0 based on w value
    T1 = 2 * pi * (m_period - 100) / w;  % Last 900 periods as effective recording period
    T0 = 2 * pi * 100 / w;               % Initialization time
    % Iterate through all D values
    for d_idx = 1:num_D
        D = D_values(d_idx);
        fprintf('\n=== Current D value: %.3f ===\n', D);
    
        QQ_values = zeros(1, num_runs);      % Store Q values for each run
        firing_rates = zeros(1, num_runs);   % Store firing rates for each run
        E_values = zeros(1, num_runs);       % Store Q values for each run

        for run = 1:num_runs
            rng(seeds(run));  % Set random seed for current run

            fprintf('Run %d/%d (Seed: %d)\n', run, num_runs, seeds(run));
        
            % Initialize variables
            v = 10; vs = 10;
            m = 0; n = 0; h = 0;

            v0 = 12 * rand();
            vs0 = rand();
            m0 = rand();
            h0 = rand();
            n0 = rand();

            firingcishu = 0;
            flag1 = 0;
            max_val = 1;
            R = 0; S = 0; Q = 0;
            E1 = 0;   % Record axonal energy
        
            %% Time loop
            Time = (0:step:T0 + T1);
            for num = 1:length(Time)
                time = Time(num);
            
                % Generate Gaussian white noise
                a = rand();
                b = rand();
                Gamma = sqrt(-4 * D * step * log(a)) .* cos(2 * pi * b);
            
                % Axonal compartment calculation
                gna = 45; Ena = 55;
                gk = 18; Ek = -80;
                gl = 0.1; El = -65;
                I_Na = -gna .* (a_m(v0, T) ./ (a_m(v0, T) + b_m(v0, T))).^3 .* h0 .* (v0 - Ena);   % Sodium current
                I_K = -gk .* m0.^4 .* (v0 - Ek);                    % Potassium current
                I_L = -gl .* (v0 - El);                           % Leak current
                I_DS = -gc .* (v0 - vs0) ./ (1 - p_);               % I_DS

                v = v0 + (I_DS + I_Na + I_K + I_L) .* step;             % Note the positions of (1-p_) and p_

                E_Na = (-I_Na) .* (v0 - Ena);                           % Sodium current energy
                E_K = (-I_K) .* (v0 - Ek);                                % Potassium current energy
                E_L = (-I_L) .* (v0 - El);                               % Leak current energy
                E_DS = (-I_DS) .* v0;                                   % I_DS energy
                E_sum = E_Na + E_K + E_L + E_DS;                            % Total synaptic energy
            
                % Dendritic compartment membrane potential calculation
                vs = vs0 + ((gc .* (v0 - vs0) / p_) - 0.1 .* (vs0 + 65) + A .* sin(w .* time)) .* step + Gamma;
            
                % Gating channel calculation
                m = m0 + (a_n(v0, T) .* (1 - m0) - b_n(v0, T) .* m0) .* step;
                n = n0 + (a_m(v0, T) .* (1 - n0) - b_m(v0, T) .* n0) .* step;
                h = h0 + (a_h(v0, T) .* (1 - h0) - b_h(v0, T) .* h0) .* step;

                % Calculate Q-value related quantities
                if time > T0 && time < T0 + T1            % Record during last 900 periods
                    E1 = E1 + E_sum;                          % Accumulate axonal energy at each time step
                    if v > -70
                        R = R + 2.0 .* v .* sin(w .* time) .* step;
                        S = S + 2.0 .* v .* cos(w .* time) .* step;
                    end
                end

                % Update variables
                v0 = v;
                vs0 = vs;
                m0 = m;
                n0 = n;
                h0 = h;

                % Calculate firing count
                if num > 20000 && v > 0.5 && flag1 == 0
                    flag1 = 1;
                end
                if flag1 == 1 && v > max_val
                    max_val = v;
                end
                if flag1 == 1 && v < max_val
                    firingcishu = firingcishu + 1;
                    flag1 = 2;
                    max_val = 1.0;
                end
                if v < -0.5 && flag1 == 2
                    flag1 = 0;
                end
            end
        
            % Calculate firing rate and Q value
            firing_rate = firingcishu / (T0 + T1 - 200) * 1000;
            firing_rates(run) = firing_rate;
        
            if firing_rate >= 0.01
                R = R / T1;
                S = S / T1;
                Q = sqrt(R^2 + S^2);
                QQ = Q / TT;
            else
                QQ = 0;
            end
            QQ_values(run) = QQ;
        
            E_values(run) = E1 ./ T1;  % Average axonal energy per unit time
        end
    
        % Calculate and store average values
        avg_Qs = mean(QQ_values);
        avg_frs = mean(firing_rates);
        avg_Es = mean(E_values);

        % Calculate efficiency
        efficient = avg_Qs ./ avg_Es;

        fprintf(fp1, '%f %f %f %f %f %f\n', w, D, avg_Qs, avg_frs, avg_Es, efficient);
    end
end
fclose(fp1);

% %% Plot Q value vs D =========================
% figure;
% plot(D_values, avg_Qs, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% xlabel('Noise Intensity D');
% ylabel('Average Q Value');
% title('Average Q Value vs Noise Intensity D');
% grid on;
% 
% %% Plot firing rate vs D
% figure;
% plot(D_values, avg_frs, 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
% xlabel('Noise Intensity D');
% ylabel('Average Firing Rate (Hz)');
% title('Average Firing Rate vs Noise Intensity D');
% grid on;
% 
% %% Plot energy vs D
% figure;
% plot(D_values, avg_Es, 'k-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% xlabel('Noise Intensity D');
% ylabel('Average E Value');
% title('Average E Value vs Noise Intensity D');
% grid on;
% 
% %% Plot energy efficiency vs D
% figure;
% plot(D_values, avg_Qs./avg_Es, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');

%% Neuron gating channel functions
function result = a_n(v, T)
    F = 2.3^((T - 23) / 10);
    result = -0.01 .* (v + 34) ./ (exp(-0.1 .* (v + 34)) - 1) .* F;
end

function result = b_n(v, T)
    F = 2.3^((T - 23) / 10);
    result = 0.125 .* exp(-(v + 44) ./ 25) .* F;
end

function result = a_m(v, T)
    F = 2.3^((T - 23) / 10);
    result = -0.1 .* (v + 33) ./ (exp(-0.1 .* (v + 33)) - 1) .* F;
end

function result = b_m(v, T)
    F = 2.3^((T - 23) / 10);
    result = 4 .* exp(-(v + 58) ./ 12) .* F;
end

function result = a_h(v, T)
    F = 2.3^((T - 23) / 10);
    result = 0.07 .* exp(-(v + 50) ./ 10) .* F;
end

function result = b_h(v, T)
    F = 2.3^((T - 23) / 10);
    result = 1 ./ (exp(-0.1 .* (v + 20)) + 1) .* F;
end