clc;clear;close all;

tic;                        %记录运行时长

M=1;                        
cm=1.0;
pi=3.1415926;

step=0.01;                     
TT=80;                        
Tao=2.0;

A=0.5;                         
w=0.01;                         

m_period=300;                   
T1=2*pi*m_period/w;             
T0=2*pi*20/w;                   

gc=0.2;                         
p_=0.5;                         
gg_=0.04;                      

fp1=fopen('.\data.txt','w');

%% D值范围和设置
D_values = 2.4:0.2:12;       
num_D = length(D_values);
avg_Qs = zeros(1, num_D);  
avg_frs = zeros(1, num_D); 
avg_Es = zeros(1, num_D);  


%% 主循环：遍历所有D值
for d_idx = 1:num_D
    D = D_values(d_idx);
    fprintf('\n=== 当前D值: %.3f ===\n', D);

    seed = 1;  
    rng(seed);  

    v=10*ones(M,TT);vs=10*ones(M,TT);

    m=zeros(M,TT); n=zeros(M,TT); h=zeros(M,TT);


    %% 生成连接矩阵============================
    K=2;
    Per=0.4;
    [matrix,x,y] = func_WS_network(TT,K,Per);


    %% 初始参数选定=================================
    I_syn=zeros(M,TT);              

    firingtime=zeros(M,TT);        
    firingcishu = 0;                
    flag1=zeros(M,TT);              
    max=ones(M,TT);                 
    E=zeros(M,TT);                  
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

        I_signal(1:TT)=A.*sin(w.*time);         

        %生成高斯白噪声
        a=rand(M,TT);                          
        b=rand(M,TT);                           
        Gamma=sqrt(-4 * D * step * log(a)) .* cos(2 * pi * b);   

        %突触部分的计算========================================================================
        gna=45; Ena=55;
        gk=18; Ek=-80;
        gl=0.1;El=-65;
        I_Na= gna.* (a_m(v0) ./ (a_m(v0) + b_m(v0))).^3 .* h0 .* (v0 - Ena) ;   
        I_K= gk.* m0.^4.* (v0-Ek);                    
        I_L= gl.* (v0 -El);                          
        I_DS= gc.* (v0 - vs0)./ (1-p_ );             

        v = v0+ (-I_DS-I_Na-I_K-I_L).* step;          

        E_Na=I_Na.*(v0 - Ena);                           
        E_K=I_K.*(v0-Ek);                               
        E_L=I_L.*(v0 -El);                               
        E_DS=I_DS.*v0;                                  
        E_sum=E_Na+E_K+E_L+E_DS;                         

        %树突部分的计算=====================================================================
        g_dl=0.1; E_dl=-65;
        I_SD= gc.*(vs0-v0)/p_;                     
        I_DL= g_dl.*(vs0-E_dl);                     

        vs= vs0+(-I_SD-I_DL+I_signal+I_syn).*step+Gamma;   

        E_ext=-(I_signal+I_syn+Gamma).*vs0;
        E_DL=I_DL.*(vs0-E_dl);                            
        E_SD=I_SD.*vs0;                                  
        ES_sum=E_DL+E_SD+E_ext;                          

        m = m0+ (a_n(v0).*(1-m0)-b_n(v0).*m0).* step;
        n = n0+ (a_m(v0).*(1-n0)-b_m(v0).*n0).* step;
        h = h0+ (a_h(v0).*(1-h0)-b_h(v0).*h0).* step;  

        if time>T0 && time<T0+T1
            E=E+(E_sum+ES_sum);                   
            if v > -72
            R =R+2.0.*v.*sin(w.*time).*step;      
            S =S+2.0.*v.*cos(w.*time).*step;      
            end
        end

        v0 = v;
        vs0 = vs;
        m0  = m;
        n0  = n;
        h0  = h;      

        flag1(v>0.5 & flag1==0)=1;
        max(flag1==1 & v>max)=v(flag1==1 & v>max);
        local=find(flag1==1 & v<max);            
        firingtime(local)=time;
        firingcishu = firingcishu + length(local);
        flag1(local)=2;
        max(local)=1;
        flag1(v<-0.5 & flag1==2)=0;        
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
        QQ=sum(Q);  
        QQ=QQ./TT;   
    else
        QQ=0;
    end
    
    % 计算并存储平均值
    avg_Qs(d_idx) = QQ;
    avg_frs(d_idx) = firing_rate;
    avg_Es(d_idx) = E_average;

    efficienct=QQ/E_average;   

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