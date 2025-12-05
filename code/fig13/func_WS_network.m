function [matrix, x, y] = func_WS_network(Num, K, Per)
% Parameter K represents the initial nearest-neighbor connections per node (i.e., the neighborhood range of the regular network)
% K is an even number, meaning during initial regular network construction, each node connects bidirectionally to K/2 nearest neighbors on each side
% Larger K indicates denser initial connections
% Parameter Per represents the random rewiring probability
% Construct small-world network
angle = 2*pi/Num:2*pi/Num:2*pi;   % Generate a circular regular network (nodes evenly distributed on the circumference)
x     = 100*sin(angle);
y     = 100*cos(angle);
matrix = zeros(Num);

% Nearest-neighbor coupling network
for i1 = 1:Num
    for j1 = i1+1:i1+K    % Each node connects to K nearest subsequent nodes
        j2 = j1;
        if j1 > Num
           j2 = mod(j1, Num);    % Modulo operation ensures indexing wraps around in the circular structure (e.g., node Num+1 becomes node 1)
        end
      matrix(i1, j2) = 1; 
      matrix(j2, i1) = 1;        % Construct bidirectional connection
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
           a = find(matrix(i1, :) == 0);   % Find positions with value 0 in row i1, locate unconnected nodes
           rand_data = randi([1, length(a)], 1, 1); 
           jjj = a(rand_data);            % Randomly select a new node
           matrix(i1, jjj) = 1;           % Establish new connection
           matrix(jjj, i1) = 1;
           matrix(i1, i1) = 0;            % Restore diagonal
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