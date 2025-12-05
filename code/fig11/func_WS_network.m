function [matrix, x, y] = func_WS_network(Num, K, Per)
% Parameter K represents the initial nearest-neighbor connections per node (i.e., the neighborhood range of the regular network)
% K is an even number, meaning each node initially connects to K/2 nearest neighbors on each side bidirectionally
% Larger K indicates denser initial connections
% Parameter Per represents the random rewiring probability
% Construct small-world network
angle = 2*pi/Num:2*pi/Num:2*pi;   % Generate a ring-shaped regular network (nodes evenly distributed on a circle)
x     = 100 * sin(angle);
y     = 100 * cos(angle);
matrix = zeros(Num);

% Nearest-neighbor coupling network
for i1 = 1:Num
    for j1 = i1+1:i1+K    % Each node connects to its K nearest subsequent neighbors
        j2 = j1;
        if j1 > Num
           j2 = mod(j1, Num);    % Modulo operation ensures circular indexing (e.g., node Num+1 wraps to node 1)
        end
        matrix(i1, j2) = 1; 
        matrix(j2, i1) = 1;        % Create bidirectional connection
    end
end

% Random rewiring process
for i1 = 1:Num
    for j1 = i1+1:i1+K
        j2 = j1;
        if j1 > Num
           j2 = mod(j1, Num);
        end
        p1 = rand();
        % Determine whether to rewire based on random probability
        if p1 < Per             % Rewire with probability Per
           matrix(i1, j2) = 0;   % Disconnect existing edge
           matrix(j2, i1) = 0;   % Disconnect edge  
           matrix(i1, i1) = inf; % inf represents infinity, temporary marker (avoid self-loop)
           a = find(matrix(i1, :) == 0);   % Find positions with value 0 in row i1 (find unconnected nodes)
           rand_data = randi([1, length(a)], 1, 1); 
           jjj = a(rand_data);            % Randomly select a new node
           matrix(i1, jjj) = 1;           % Establish new connection
           matrix(jjj, i1) = 1;
           matrix(i1, i1) = 0;            % Restore diagonal element
        end
    end
end

% % Visualize network
% figure;
% G = graph(matrix);
% p = plot(G, 'XData', x, 'YData', y, 'NodeColor', 'b', 'EdgeColor', [0.5 0.5 0.5], 'MarkerSize', 6);
% title(sprintf('Small-world network (N=%d, K=%d, p=%.2f)', Num, K, Per));
% axis equal;
% end