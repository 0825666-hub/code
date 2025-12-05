clc; clear; close all;

TT = 1; % Single neuron
cm = 1.0;
pi = 3.1415926;

step = 0.01; % Time step
Tao = 2;

v = 0; vs = 0;
m = 0; n = 0; h = 0;

fp1 = fopen('.\Q.dat', 'w');
fp2 = fopen('.\SynchronizationFactor.dat', 'w');

%% Loop for w===================================================
for w_num = 1:46
    fprintf('%d\n', w_num); % Display progress
    w = 0.01 + 0.002 * (w_num - 1); % w ranges from 0.01 to 0.10

    m_period = 20;
    T1 = 2 * pi * m_period / w;
    T0 = 200;
    D = 0; % Diffusion coefficient, i.e., noise intensity (adjustable as needed)
    T = 23; % Temperature factor temporarily not added

    %% Loop for A===========================================================
    for A_num = 1:85
        fprintf('%d\n', A_num); % Display progress
        A = 0.015 * (A_num - 1); % A ranges from 0 to 1.26

        time = 0;

        gc = 0.5; % Conductance between two compartments
        p_ = 0.7; % Ratio of dendritic area to total cell area
        firingrate = 0;
        firingcishu = 0;
        v0 = 0; vs0 = 0;
        m0 = 0; h0 = 0; n0 = 0;
        flag1 = 0; % Variable to assist in recording spikes
        firingtime = 0;
        max_val = 1; % Variable to assist in recording spikes
        R = 0; S = 0; Q = 0; % Values needed to calculate Q

        H = 0;
        bb_ = 0;

        ran2 = rand(); ran3 = rand(); ran4 = rand(); % Generate three random numbers
        v0 = 12 * ran2;
        vs0 = ran2;
        m0 = ran2;
        h0 = ran3;
        n0 = ran4; % Parameter initialization

        %% Time loop=================================================
        Time = (0:step:T0+T1);
        for num = 1:length(Time)
            time = Time(num);

            % Generate Gaussian white noise
            a = rand(); % Generate random number a
            b = rand(); % Generate random number b
            Gamma = sqrt(-4 * D * step * log(a)) .* cos(2 * pi * b); % Calculate Gamma

            v = v0 + (-gc .* (v0 - vs0) ./ (1 - p_) ...
                - 45 .* (a_m(v0, T) ./ (a_m(v0, T) + b_m(v0, T))).^3 .* h0 .* (v0 - 55) ...
                - 18 .* m0.^4 .* (v0 + 80) - 0.1 .* (v0 + 65)) .* step; 
            vs = vs0 + ((gc .* (v0 - vs0) / p_) - 0.1 .* (vs0 + 65) + A .* sin(w .* time) + Gamma) .* step; 
            m = m0 + (a_n(v0, T) .* (1 - m0) - b_n(v0, T) .* m0) .* step;
            n = n0 + (a_m(v0, T) .* (1 - n0) - b_m(v0, T) .* n0) .* step;
            h = h0 + (a_h(v0, T) .* (1 - h0) - b_h(v0, T) .* h0) .* step;

            if time > T0 && time < T0 + T1
                R = R + 2.0 .* v .* sin(w .* time) .* step;
                S = S + 2.0 .* v .* cos(w .* time) .* step;
            end

            v0 = v;
            vs0 = vs;
            m0 = m;
            n0 = n;
            h0 = h;

            if num > 20000 && v > 0.5 && flag1 == 0
                flag1 = 1;
            end
            if flag1 == 1 && v > max_val
                max_val = v;
            end
            if flag1 == 1 && v < max_val
                firingcishu = firingcishu + 1;
                flag1 = 2;
                firingtime = time;
                max_val = 1.0;
                bb_ = bb_ + 1;
            end
            if v < -0.5
                if flag1 == 2
                    flag1 = 0;
                end
            end
        end

        firingrate = firingcishu / (T0 + T1 - 200) * 1000;
        %% Calculation of Q=============================================================

        QQ = 0;

        R = R ./ T1;
        S = S ./ T1;
        Q = sqrt(R.^2 + S.^2);
        QQ = Q;

        QQ = QQ ./ TT;

        H = H ./ num;

        fprintf(fp1, '%f %f %f %f \n', D, QQ, H, H ./ QQ);
        fprintf(fp2, '%f %f %f \n', w, A, bb_);
    end
end

fclose(fp1);
fclose(fp2);

%% Data processing ======================================

fp_data = load('SynchronizationFactor.dat');
fp_bb = fp_data(:, 3);
Z = reshape(fp_bb, length(0:0.015:1.26), length(0.01:0.002:0.1));

figure;
contourf(0.01:0.002:0.1, 0:0.015:1.26, Z, 10); % Filled contour plot, 10 indicates number of contour lines
xlabel('w'); ylabel('A');
title('Firing Range');
xlim([0.01 0.08]);
ylim([0.6 1.26]);
colorbar;

%% Firing Area Plot===========================================
% Categorize the matrix
firing_area = Z ~= 0; % Non-zero area is 1, zero area is 0

% Plot heatmap
imagesc(0.01:0.002:0.1, 0:0.015:1.26, firing_area);
colormap([0.4 0.6 0.8; 1 0.4 0.4]); % Color differentiation
xlabel('w'); ylabel('A');
title('Firing Range');
xlim([0.01 0.085]);
ylim([0.6 1.26]);

% Set Y-axis to normal direction (bottom to top)
set(gca, 'YDir', 'normal'); % Y-axis from bottom to top

% Proportion of firing area======================
% Count the number of zero elements in the array
num_zeros = sum(fp_bb == 0);

% Calculate total number of elements in the array
total_elements = numel(fp_bb);

% Calculate the proportion of zeros
zero_ratio = num_zeros / total_elements;

title(['Firing Range Proportion: ', num2str(1 - zero_ratio)]);

%% ================================================
function result = a_n(v, T)
    F = 2.3^((T - 23) / 10); % Add temperature factor F
    result = -0.01 .* (v + 34) ./ (exp(-0.1 .* (v + 34)) - 1) .* F;
end

function result = b_n(v, T)
    F = 2.3^((T - 23) / 10); % Add temperature factor F
    result = 0.125 .* exp(-(v + 44) ./ 25) .* F;
end

function result = a_m(v, T)
    F = 2.3^((T - 23) / 10); % Add temperature factor F
    result = -0.1 .* (v + 33) ./ (exp(-0.1 .* (v + 33)) - 1) .* F;
end

function result = b_m(v, T)
    F = 2.3^((T - 23) / 10); % Add temperature factor F
    result = 4 .* exp(-(v + 58) ./ 12) .* F;
end

function result = a_h(v, T)
    F = 2.3^((T - 23) / 10); % Add temperature factor F
    result = 0.07 .* exp(-(v + 50) ./ 10) .* F;
end

function result = b_h(v, T)
    F = 2.3^((T - 23) / 10); % Add temperature factor F
    result = 1 ./ (exp(-0.1 .* (v + 20)) + 1) .* F;
end