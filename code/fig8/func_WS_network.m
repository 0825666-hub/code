function [matrix, x, y] = func_WS_network(Num, K, Per)

%rng(1);
% Construct small-world network
angle = 2*pi/Num:2*pi/Num:2*pi;  
x     = 100 * sin(angle);
y     = 100 * cos(angle);
matrix = zeros(Num);

% Nearest-neighbor coupling network
for i1 = 1:Num
    for j1 = i1+1:i1+K    
        j2 = j1;
        if j1 > Num
           j2 = mod(j1, Num);    
        end
        matrix(i1, j2) = 1; 
        matrix(j2, i1) = 1;        
    end
end

% Rewiring edges with probability Per
for i1 = 1:Num
    for j1 = i1+1:i1+K
        j2 = j1;
        if j1 > Num
           j2 = mod(j1, Num);
        end
        p1 = rand();
        
        if p1 < Per             
           % Disconnect existing edge
           matrix(i1, j2) = 0;   
           matrix(j2, i1) = 0;   
           
           % Mark self-connection as infinite to avoid reconnecting to itself
           matrix(i1, i1) = inf; 
           
           % Find all nodes not currently connected to i1
           a             = find(matrix(i1, :) == 0);  
           
           % Randomly select one node to connect to
           rand_data     = randi([1, length(a)], 1, 1); 
           jjj           = a(rand_data);            
           
           % Create new connection
           matrix(i1, jjj) = 1;                       
           matrix(jjj, i1) = 1;
           
           % Reset self-connection
           matrix(i1, i1) = 0;                       
        end
    end
end

% % Visualize the network
% figure;
% G = graph(matrix);
% p = plot(G, 'XData', x, 'YData', y, 'NodeColor', 'b', 'EdgeColor', [0.5 0.5 0.5], 'MarkerSize', 6);
% title(sprintf('Small-world network (N=%d, K=%d, p=%.2f)', Num, K, Per));
% axis equal;
end
