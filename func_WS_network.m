function [matrix,x,y] = func_WS_network(Num,K,Per)
%参数 K 表示 每个节点的初始最近邻连接数（即规则网络的邻域范围）
%K 是偶数，代表在构建初始规则网络时，每个节点会与左右各 K/2 个最近邻节点双向连接
%K 越大，初始连接越密集
%Per 表示 随机重连概率
rng(1);
% 构建小世界 
angle = 2*pi/Num:2*pi/Num:2*pi;   %生成一个环形规则网络（节点均匀分布在圆周上）
x     = 100*sin(angle);
y     = 100*cos(angle);
matrix= zeros(Num);

% 最近邻耦合网络
for i1=1:Num
    for j1=i1+1:i1+K    %每个节点连接其后的 K 个最近邻节点
        j2=j1;
        if j1 > Num
           j2 = mod(j1,Num);    % mod取余运算，确保索引在环形结构中循环（如第 Num+1 个节点即第 1 个节点）
        end
      matrix(i1,j2) = 1; 
      matrix(j2,i1) = 1;        %构建双向连接
    end
end

for i1=1:Num
    for j1=i1+1:i1+K
        j2=j1;
        if j1>Num
           j2=mod(j1,Num);
        end
        p1 = rand();
        % 根据随机概率判断是否进行连边
        if p1 < Per             % 以概率Per重连
           matrix(i1,j2) = 0;   % 断开原有边
           matrix(j2,i1) = 0;   % 断边  
           matrix(i1,i1) = inf; % inf表示无穷大，临时标记（避免自环）
           a             = find(matrix(i1,:)==0);   % 找到i1行中为0的位置，找到未连接的节点
           rand_data     = randi([1,length(a)],1,1); 
           jjj           = a(rand_data);            % 随机选择一个新节点
           matrix(i1,jjj)= 1;                       % 建立新连接
           matrix(jjj,i1)= 1;
           matrix(i1,i1) = 0;                       % 恢复对角线
        end
    end
end

% % 可视化网络
% figure;
% G = graph(matrix);
% p = plot(G, 'XData', x, 'YData', y, 'NodeColor', 'b', 'EdgeColor', [0.5 0.5 0.5], 'MarkerSize', 6);
% title(sprintf('小世界网络 (N=%d, K=%d, p=%.2f)', Num, K, Per));
% axis equal;
% end