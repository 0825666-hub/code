clc; clear; close all;

%% Parameter Settings
TT = 1;
cm = 1.0;
pi = 3.1415926;

step = 0.01;
Tao = 2;

w = 0.01;
A = 0.5;   % A cannot be too large       
T = 23;    % Temperature factor temporarily not added      

m_period = 100;   
T1 = 2 * pi * m_period / w;
T0 = 2 * pi * 10 / w;

gc = 0.2;         % Conductance between two compartments
p_ = 0.7;         % Ratio of dendritic area to total cell area

%% D value range and settings
D_values = 0:0.25:20;  
num_D = length(D_values);
avg_Qs = zeros(1, num_D);  
avg_frs = zeros(1, num_D); 
avg_Es = zeros(1, num_D);  

%% Run multiple simulations
num_runs = 10;  
seeds = [1,2,3,4,5,6,7,8,9,10];  

%% Main loop: iterate through all D values
for d_idx = 1:num_D
    D = D_values(d_idx);
    fprintf('\n=== Current D value: %.3f ===\n', D);
    
    QQ_values = zeros(1, num_runs);      
    firing_rates = zeros(1, num_runs);   
    E_values = zeros(1, num_runs);      

    for run = 1:num_runs
        rng(seeds(run));  
        
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
        E1 = 0;   
        
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
            I_Na = -gna .* (a_m(v0, T) ./ (a_m(v0, T) + b_m(v0, T))).^3 .* h0 .* (v0 - Ena);   
            I_K = -gk .* m0.^4 .* (v0 - Ek);                    
            I_L = -gl .* (v0 - El);                          
            I_DS = -gc .* (v0 - vs0) ./ (1 - p_);               

            v = v0 + (I_DS + I_Na + I_K + I_L) .* step;             

            E_Na = (-I_Na) .* (v0 - Ena);                           
            E_K = (-I_K) .* (v0 - Ek);                                
            E_L = (-I_L) .* (v0 - El);                               
            E_DS = (-I_DS) .* v0;                                   
            E_sum = E_Na + E_K + E_L + E_DS;                            
            
            % Dendritic compartment membrane potential calculation
            vs = vs0 + ((gc .* (v0 - vs0) / p_) - 0.1 .* (vs0 + 65) + A .* sin(w .* time)) .* step + Gamma;
            
            % Gating channel calculation
            m = m0 + (a_n(v0, T) .* (1 - m0) - b_n(v0, T) .* m0) .* step;
            n = n0 + (a_m(v0, T) .* (1 - n0) - b_m(v0, T) .* n0) .* step;
            h = h0 + (a_h(v0, T) .* (1 - h0) - b_h(v0, T) .* h0) .* step;
            
            % Calculate Q-value related quantities
            if time > T0 && time < T0 + T1 
                E1 = E1 + E_sum;         
                if v > -70
                    R = R + 2.0 .* v .* sin(w .* time) .* step;
                    S = S + 2.0 .* v .* cos(w .* time) .* step;
                end  
            end
            
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
        
        % Calculate firing rate and Q-value
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
        
        E_values(run) = E1 ./ T1;  
    end
    
    % Calculate average values for this D
    avg_Qs(d_idx) = mean(QQ_values);
    avg_frs(d_idx) = mean(firing_rates);
    avg_Es(d_idx) = mean(E_values);
    
    fprintf('Average Q value: %.4f\n', avg_Qs(d_idx));
    fprintf('Average firing rate: %.4f Hz\n', avg_frs(d_idx));
    fprintf('Average energy: %.4f\n', avg_Es(d_idx));
end

%% ===========================================
% Find trough locations to set initial values to zero
[~, trough_locs1] = findpeaks(-avg_Qs);
[~, trough_locs2] = findpeaks(-avg_frs);
[~, trough_locs3] = findpeaks(-avg_Es);

avg_Qs(1:trough_locs1-1) = 0;
avg_frs(1:trough_locs2-1) = 0;
avg_Es(1:trough_locs3-1) = 0;

% Plot Q value vs D =========================
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