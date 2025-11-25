clc;clear;close all;

tic;                        %记录运行时长

M=1;                        % 层数
cm=1.0;
pi=3.1415926;

step=0.01;                     %时间步长
TT=80;                         %每层神经元的数目
Tao=2.0;

A=0.5;                          %外部弱信号的强度
w=0.01;                         %外部弱信号的频率

m_period=300;                   %模拟运行的周期个数
T1=2*pi*m_period/w;             %记录期间的运行时长
T0=2*pi*20/w;                   %初始瞬态（无需记录的时间）

gc=0.2;                         %双室之间的突触电导
p_=0.5;                         %树突面积比
gg_=0.04;                       %节点之间的电导（耦合强度）

fp1=fopen('.\data.txt','w');

%% D值范围和设置
D_values = 2.4:0.2:12;       % D从0到12，步长0.02,60个数据点
num_D = length(D_values);
avg_Qs = zeros(1, num_D);  % 存储每个D对应的平均Q值
avg_frs = zeros(1, num_D); % 存储每个D对应的平均放电率
avg_Es = zeros(1, num_D);  % 存储每个D对应的平均能量E


%% 主循环：遍历所有D值
for d_idx = 1:num_D
    D = D_values(d_idx);
    fprintf('\n=== 当前D值: %.3f ===\n', D);

    seed = 1;  % 1个随机种子
    rng(seed);  % 设置当前运行的随机种子

    v=10*ones(M,TT);vs=10*ones(M,TT);

    m=zeros(M,TT); n=zeros(M,TT); h=zeros(M,TT);


    %% 生成连接矩阵============================
    K=2;
    Per=0.4;
    [matrix,x,y] = func_WS_network(TT,K,Per);


    %% 初始参数选定=================================
    I_syn=zeros(M,TT);              %节间电流

    firingtime=zeros(M,TT);         %记录放电时间
    firingcishu = 0;                %记录总放电数量的量
    flag1=zeros(M,TT);              %辅助判断放电的量
    max=ones(M,TT);                 %辅助判断放电的量
    E=zeros(M,TT);                  %能量
    R=zeros(M,TT);   S=zeros(M,TT);   Q=zeros(M,TT);

    v0=12*rand(M,TT);
    vs0=rand(M,TT);
    m0=rand(M,TT);
    h0=rand(M,TT);
    n0=rand(M,TT);

    %% 时间循环=============================================
    Time=(0:step:T0+T1);
    for num=1:length(Time)
        time=Time(num);
        if mod(num,100000)==0
            fprintf('%d\n',num);
        end

        if time>T0
            alpha_fun=alpha(time-firingtime);
            IsynE=-matrix.*gg_.*alpha_fun.*(v-0);
            I_syn=sum(IsynE,2).';                %转置
        end

        I_signal(1:TT)=A.*sin(w.*time);          %外部输入信号

        %生成高斯白噪声
        a=rand(M,TT);                           %生成随机数 a
        b=rand(M,TT);                           %生成随机数 b
        Gamma=sqrt(-4 * D * step * log(a)) .* cos(2 * pi * b);   % 计算Gamma

        %突触部分的计算========================================================================
        gna=45; Ena=55;
        gk=18; Ek=-80;
        gl=0.1;El=-65;
        I_Na= gna.* (a_m(v0) ./ (a_m(v0) + b_m(v0))).^3 .* h0 .* (v0 - Ena) ;   %钠电流
        I_K= gk.* m0.^4.* (v0-Ek);                    %钾电流
        I_L= gl.* (v0 -El);                           %漏电流
        I_DS= gc.* (v0 - vs0)./ (1-p_ );              %I_DS

        v = v0+ (-I_DS-I_Na-I_K-I_L).* step;          %注意（1-p_）和p_的位置

        E_Na=I_Na.*(v0 - Ena);                           %钠电流能量
        E_K=I_K.*(v0-Ek);                                %钾电流能量
        E_L=I_L.*(v0 -El);                               %漏电流能量
        E_DS=I_DS.*v0;                                   %I_DS能量
        E_sum=E_Na+E_K+E_L+E_DS;                         %突触总能量

        %树突部分的计算=====================================================================
        g_dl=0.1; E_dl=-65;
        I_SD= gc.*(vs0-v0)/p_;                      %I_SD
        I_DL= g_dl.*(vs0-E_dl);                     %漏电流

        vs= vs0+(-I_SD-I_DL+I_signal+I_syn).*step+Gamma;   

        E_ext=-(I_signal+I_syn+Gamma).*vs0;
        E_DL=I_DL.*(vs0-E_dl);                            %漏电流能量
        E_SD=I_SD.*vs0;                                   %I_SD能量
        ES_sum=E_DL+E_SD+E_ext;                           %树突总能量

        m = m0+ (a_n(v0).*(1-m0)-b_n(v0).*m0).* step;
        n = n0+ (a_m(v0).*(1-n0)-b_m(v0).*n0).* step;
        h = h0+ (a_h(v0).*(1-h0)-b_h(v0).*h0).* step;  %门控变量

        if time>T0 && time<T0+T1
            E=E+(E_sum+ES_sum);                   %能量
            if v > -72
            R =R+2.0.*v.*sin(w.*time).*step;      %Qsin
            S =S+2.0.*v.*cos(w.*time).*step;      %Qcos
            end
        end

        v0 = v;
        vs0 = vs;
        m0  = m;
        n0  = n;
        h0  = h;      %数据迭代

        flag1(v>0.5 & flag1==0)=1;
        max(flag1==1 & v>max)=v(flag1==1 & v>max);
        local=find(flag1==1 & v<max);            %放电的节点
        firingtime(local)=time;
        firingcishu = firingcishu + length(local);
        flag1(local)=2;
        max(local)=1;
        flag1(v<-0.5 & flag1==2)=0;        %reset重置
    end

    %% 计算Q和R======================================================
    % 计算放电率
    firing_rate = (firingcishu/TT)/ (T0 + T1 - 200) * 1000;

    %计算能量
    E_average=(sum(E)/TT)/T1;

    if firing_rate >= 0.01
        R=R./T1;
        S=S./T1;
        Q=sqrt(R.^2 + S.^2);
        QQ=sum(Q);   %求和
        QQ=QQ./TT;   %求平均值
    else
        QQ=0;
    end
    
    % 计算并存储平均值
    avg_Qs(d_idx) = QQ;
    avg_frs(d_idx) = firing_rate;
    avg_Es(d_idx) = E_average;

    efficienct=QQ/E_average;   %能量效率

    fprintf(fp1,'%f %f %f %f %f\n',D,avg_Qs(d_idx),avg_frs(d_idx),avg_Es(d_idx));
end

fclose(fp1);

runtime=toc;  %记录运行时长

%% =======================================================
clc;clear;close all;
fp1=fopen('.\data1.txt','w');
data=load('data.txt');
col1=data(:,1);
col2=data(:,2)*0.4;
col3=data(:,3);
col4=data(:,4)*2.2;
col5=col2./col4;
for i=1:length(col1)
    fprintf(fp1,'%f %f %f %f %f\n',col1(i),col2(i),col3(i),col4(i),col5(i));
end
fclose(fp1);

%% 所用到的函数================================================
function result = a_n(v)
    result = -0.01 .* (v + 34) ./ (exp(-0.1 .* (v + 34)) - 1);
end

function result = b_n(v)
    result = 0.125 .* exp(-(v + 44) ./ 25);
end

function result = a_m(v)
    result = -0.1 .* (v + 33) ./ (exp(-0.1 .* (v + 33)) - 1);
end

function result = b_m(v)
    result = 4 .* exp(-(v + 58) ./ 12);
end

function result = a_h(v)
    result = 0.07 .* exp(-(v + 50) ./ 10);
end

function result = b_h(v)
    result = 1 ./ (exp(-0.1 .* (v + 20)) + 1);
end

function result = alpha(t)
    Tao=2.0;
    result = t .* exp(-t ./ Tao) ./ Tao;
end