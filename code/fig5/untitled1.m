clc; clear; close all;

TT = 1;           % Number of neurons (single neuron)
cm = 1.0;
pi = 3.1415926;

step = 0.01; % Time step size
Tao = 2;

v = 0; vs = 0;
m = 0; n = 0; h = 0;

fp1 = fopen('.\v.dat','w');

w = 0.05;   % If w is too large, it will cause the neuron refractory period
A = 0.6;
T = 23;     % Temperature factor is temporarily not added

m_period = 100;   % Value is related to w
T1 = 2 * pi * m_period / w;
T0 = 200;         % Initialization time before recording

time = 0;

gc = 0.1;         % Conductance between the two compartments
p_ = 0.8;         % Ratio of cell area occupied by the dendrite
firingrate = 0;
firingcishu = 0;
v0 = 0; vs0 = 0;
m0 = 0; h0 = 0; n0 = 0;
flag1 = 0;     % Helper variable to record spikes
firingtime = 0;
max = 1;       % Helper variable to record spikes
R = 0; S = 0; Q = 0; % Values needed to calculate Q
E1 = 0;   % Record axon energy
E2 = 0;   % Record dendrite energy
E1_data = zeros(length(T0:step:T0+T1), 1);

H = 0;

ran2 = rand(); ran3 = rand(); ran4 = rand();   % Generate three random numbers
v0 = 12 * ran2;
vs0 = ran2;
m0 = ran2;
h0 = ran3;
n0 = ran4;                                % Parameter initialization

D = 0.8;                                    % Diffusion coefficient, i.e., noise intensity (adjustable)

%% Main time loop =================================================
Time = (0:step:T0+T1);
for num = 1:length(Time)
    time = Time(num);
    if mod(num, 100000) == 0
        fprintf('%d\n', num);
    end

    % Generate Gaussian white noise
    a = rand();                           % Generate random number a
    b = rand();                           % Generate random number b
    Gamma = sqrt(-4 * D * step * log(a)) .* cos(2 * pi * b);   % Calculate Gamma

    % Calculations for the synaptic (soma) compartment ==========================================
    gna = 45; Ena = 55;
    gk = 18; Ek = -80;
    gl = 0.1; El = -65;
    I_Na = - gna .* (a_m(v0, T) ./ (a_m(v0, T) + b_m(v0, T))).^3 .* h0 .* (v0 - Ena);   % Sodium current
    I_K = - gk .* m0.^4 .* (v0 - Ek);                    % Potassium current
    I_L = - gl .* (v0 - El);                           % Leak current
    I_DS = -gc .* (v0 - vs0) ./ (1 - p_);               % I_DS (current from dendrite to soma)

    v = v0 + (I_DS + I_Na + I_K + I_L) .* step;        % Note the positions of (1-p_) and p_

    E_Na = (-I_Na) .* (v0 - Ena);                           % Sodium current energy
    E_K = (-I_K) .* (v0 - Ek);                                % Potassium current energy
    E_L = (-I_L) .* (v0 - El);                               % Leak current energy
    E_DS = (-I_DS) .* v0;                                   % I_DS energy
    E_sum = E_Na + E_K + E_L + E_DS;                            % Total energy in the somatic compartment

    % Calculations for the dendritic compartment ===========================================
    g_dl = 0.1; E_dl = -65;
    I_SD = gc .* (v0 - vs0) / p_;                            % I_SD (current from soma to dendrite)
    I_DL = -g_dl .* (vs0 - E_dl);                          % Leak current

    vs = vs0 + (I_SD + I_DL + A .* sin(w .* time)) .* step + Gamma;

    % E_ext = (A .* sin(w .* time) + Gamma) .* vs0;               % Energy from external current (to be included?)
    % E_DL = I_DL .* (vs0 - E_dl);                            % Leak current energy
    % E_SD = I_SD .* vs0;                                   % I_SD energy
    % ES_sum = E_DL + E_SD + E_ext;                           % Total energy in the dendritic compartment

    % Gating variable calculations ==================================================================
    m = m0 + (a_n(v0, T) .* (1 - m0) - b_n(v0, T) .* m0) .* step;
    n = n0 + (a_m(v0, T) .* (1 - n0) - b_m(v0, T) .* n0) .* step;
    h = h0 + (a_h(v0, T) .* (1 - h0) - b_h(v0, T) .* h0) .* step;

    % Data recording during the measurement period
    if time > T0 && time < T0 + T1
        R = R + 2.0 .* v .* sin(w .* time) .* step;
        S = S + 2.0 .* v .* cos(w .* time) .* step;
        E1 = E1 + E_sum;         % Accumulate axon energy at each time step
        % E2 = E2 + ES_sum;        % Accumulate dendrite energy at each time step
        E1_data(num - T0 ./ step) = E_sum;   % Record axon energy at each time step
    end

    % Update state variables for the next iteration
    v0 = v;
    vs0 = vs;
    m0 = m;
    n0 = n;
    h0 = h;

    % Spike detection logic
    if num > 20000 && v > 0.5 && flag1 == 0
        flag1 = 1;
    end
    if flag1 == 1 && v > max
        max = v;
    end
    if flag1 == 1 && v < max
        firingcishu = firingcishu + 1;
        flag1 = 2;
        firingtime = time;
        max = 1.0;
    end
    if v < -0.5
        if flag1 == 2
            flag1 = 0;
        end
    end

    % Write data to file at specified intervals
    if mod(num, 30) == 0 && num > 500000 && num < 2500000
        fprintf(fp1, '%f %f %f %f\n', ...
            time, v, v, 10 * sin(w * time));
    end
end

% Calculate firing rate
firingrate = firingcishu / (T0 + T1 - 200) * 1000;

QQ = 0;

% Calculate signal amplitude Q
R = R ./ T1;
S = S ./ T1;
Q = sqrt(R.^2 + S.^2);
QQ = Q;

QQ = QQ ./ TT;

% Calculate average energy consumption
E1_average = E1 ./ T1;  % Average axon energy per unit time
E2_average = E2 ./ T1;  % Average dendrite energy per unit time

fclose(fp1);

%% Plotting section =============================================================
% Load data and plot voltage trace
fp_data = load('.\v.dat');
time = fp_data(:, 1);
v = fp_data(:, 2);
signal = fp_data(:, 4);
figure;
plot(time, v, 'k', 'LineWidth', 1);
xlim([10000 12500]);
hold on;
plot(time, signal, 'r', 'LineWidth', 2);
xlim([10000 12500]);
hold off;

%% Plotting section for energy =============================================================
figure;
plot(T0:step:T0+T1, E1_data);
xlim([10000 12500]);

%% ================================================
% Gating variable rate functions with temperature factor
function result = a_n(v, T)
    F = 2.3^((T - 23) / 10);     % Temperature scaling factor F
    result = -0.01 .* (v + 34) ./ (exp(-0.1 .* (v + 34)) - 1) .* F;
end

function result = b_n(v, T)
    F = 2.3^((T - 23) / 10);     % Temperature scaling factor F
    result = 0.125 .* exp(-(v + 44) ./ 25) .* F;
end

function result = a_m(v, T)
    F = 2.3^((T - 23) / 10);     % Temperature scaling factor F
    result = -0.1 .* (v + 33) ./ (exp(-0.1 .* (v + 33)) - 1) .* F;
end

function result = b_m(v, T)
    F = 2.3^((T - 23) / 10);     % Temperature scaling factor F
    result = 4 .* exp(-(v + 58) ./ 12) .* F;
end

function result = a_h(v, T)
    F = 2.3^((T - 23) / 10);     % Temperature scaling factor F
    result = 0.07 .* exp(-(v + 50) ./ 10) .* F;
end

function result = b_h(v, T)
    F = 2.3^((T - 23) / 10);     % Temperature scaling factor F
    result = 1 ./ (exp(-0.1 .* (v + 20)) + 1) .* F;
end
