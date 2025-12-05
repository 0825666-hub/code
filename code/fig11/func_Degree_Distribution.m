function [Dds, Dds_avg, M, P_Dds] = func_Degree_Distribution(matrix)
 
Num = size(matrix, 2); % Returns the number of columns
Dds = zeros(1, Num);   % Create a horizontal array

for i = 1:Num
    Dds(i) = sum(matrix(i, :)); % Get degree of each node
end
Dds_avg = mean(Dds);         % Calculate average degree
 
M = max(Dds);                % Maximum degree

Num_Dds = zeros(1, M+1);
for i = 1:M+1   
    Num_Dds(i) = length(find(Dds == i-1)); % Find count of each degree value
end
P_Dds = zeros(1, M+1);
P_Dds(:) = Num_Dds(:) ./ sum(Num_Dds);     % Degree distribution