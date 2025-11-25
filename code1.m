clc; clear; close all;

%% 参数设置
TT = 1;
cm = 1.0;
pi = 3.1415926;

step = 0.01;
Tao = 2;

w = 0.01;
A = 0.5;          
T = 23;          

m_period = 100;   
T1 = 2 * pi * m_period / w;
T0 = 2 * pi * 10 / w;

gc = 0.2;         
p_ = 0.7;         

%% D值范围和设置
D_values = 0:0.25:20;  
num_D = length(D_values);
avg_Qs = zeros(1, num_D);  
avg_frs = zeros(1, num_D); 
avg_Es = zeros(1, num_D);  

%% 运行多次模拟
num_runs = 10;  
seeds = [1,2,3,4,5,6,7,8,9,10];  

%% 主循环：遍历所有D值
for d_idx = 1:num_D
    D = D_values(d_idx);
    fprintf('\n=== 当前D值: %.3f ===\n', D);
    
    QQ_values = zeros(1, num_runs);      
    firing_rates = zeros(1, num_runs);   
    E_values = zeros(1, num_runs);      

    for run = 1:num_runs
        rng(seeds(run));  
        
        fprintf('Run %d/%d (Seed: %d)\n', run, num_runs, seeds(run));
        
        % 初始化变量
        v = 10; vs = 10;
        m = 0; n = 0; h = 0;
        
        v0 = 12 * rand();
        vs0 = rand();
        m0 = rand();
        h0 = rand();
        n0 = rand();
        
        firingcishu = 0;
        flag1 = 0;
        max = 1;
        R = 0; S = 0; Q = 0;
        E1=0;   
        
        %% 时间循环
        Time = (0:step:T0 + T1);
        for num = 1:length(Time)
            time = Time(num);
            
            
            a = rand();
            b = rand();
            Gamma = sqrt(-4 * D * step * log(a)) .* cos(2 * pi * b);
            
            %轴突部分的计算
            gna=45; Ena=55;
            gk=18; Ek=-80;
            gl=0.1;El=-65;
            I_Na= - gna.* (a_m(v0,T) ./ (a_m(v0,T) + b_m(v0,T))).^3 .* h0 .* (v0 - Ena) ;   
            I_K= - gk.* m0.^4.* (v0-Ek);                    
            I_L= - gl.* (v0 -El);                          
            I_DS= -gc.* (v0 - vs0)./ (1-p_ );               

            v = v0+ (I_DS+I_Na+I_K+I_L).* step;             

            E_Na=(-I_Na).*(v0 - Ena);                           
            E_K=(-I_K).*(v0-Ek);                                
            E_L=(-I_L).*(v0 -El);                               
            E_DS=(-I_DS).*v0;                                   
            E_sum=E_Na+E_K+E_L+E_DS;                            
            
            %树突部分的膜电位计算
            vs = vs0 + ((gc .* (v0 - vs0) / p_) - 0.1 .* (vs0 + 65) + A .* sin(w .* time)) .* step + Gamma;
            
            %门控通道的计算
            m = m0 + (a_n(v0, T) .* (1 - m0) - b_n(v0, T) .* m0) .* step;
            n = n0 + (a_m(v0, T) .* (1 - n0) - b_m(v0, T) .* n0) .* step;
            h = h0 + (a_h(v0, T) .* (1 - h0) - b_h(v0, T) .* h0) .* step;
            
           
            if time > T0 && time < T0 + T1 
                E1=E1+E_sum;         
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
            
            
            if num > 20000 && v > 0.5 && flag1 == 0
                flag1 = 1;
            end
            if flag1 == 1 && v > max
                max = v;
            end
            if flag1 == 1 && v < max
                firingcishu = firingcishu + 1;
                flag1 = 2;
                max = 1.0;
            end
            if v < -0.5 && flag1 == 2
                flag1 = 0;
            end
        end
        
        
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
        
        E_values(run)=E1./T1;  
    end
    
    
    avg_Qs(d_idx) = mean(QQ_values);
    avg_frs(d_idx) = mean(firing_rates);
    avg_Es(d_idx) = mean(E_values);
    
    fprintf('平均 Q 值: %.4f\n', avg_Qs(d_idx));
    fprintf('平均放电率: %.4f Hz\n', avg_frs(d_idx));
    fprintf('平均能量: %.4f\n', avg_Es(d_idx));
end

%% ===========================================
[~, trough_locs1] = findpeaks(-avg_Qs);
[~, trough_locs2] = findpeaks(-avg_frs);
[~, trough_locs3] = findpeaks(-avg_Es);

avg_Qs(1:trough_locs1-1)=0;
avg_frs(1:trough_locs2-1)=0;
avg_Es(1:trough_locs3-1)=0;

%绘制Q值随D变化的图像=========================
% figure;
% plot(D_values, avg_Qs, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% xlabel('噪声强度 D');
% ylabel('平均 Q 值');
% title('平均 Q 值随噪声强度 D 的变化');
% grid on;
% 
% %% 绘制放电率随D变化的图像
% figure;
% plot(D_values, avg_frs, 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
% xlabel('噪声强度 D');
% ylabel('平均放电率 (Hz)');
% title('平均放电率随噪声强度 D 的变化');
% grid on;
% 
% %% 绘制能量随D变化的图像
% figure;
% plot(D_values, avg_Es, 'k-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% xlabel('噪声强度 D');
% ylabel('平均 E 值');
% title('平均 E 值随噪声强度 D 的变化');
% grid on;
% 
% %% 绘制能量效率随D变化的图像
% figure;
% plot(D_values, avg_Qs./avg_Es, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
%% 神经元门控通道函数
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