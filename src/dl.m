function [U, dist_array, avg_dist_array, count_inactive] = dl(P, epsilon, stopping_func)
% DL Learns a dictionary for an input dataset using the distance from
% line-segment pairs algorithm.
%
%   INPUT:
%   P             - the input dataset as a matrix with columns as
%                   datapoints and rows as dimensions
%   epsilon       - error tolerance for each datapoint
%   stopping_func - either @max or @mean which is to be used as the
%                   stopping criterion (default = @max)
%
%   OUTPUT:
%   U              - the dictionary learned by the algorithm as a matrix
%   dist_array     - distance of farthest point from the pairs of line
%                    segments induced by dictionary atoms at each iteration
%                    of the algorithm
%   avg_dist_array - average distance points from the pairs of line
%                    segments induced by dictionary atoms at each iteration
%                    of the algorithm
%   count_inactive - number of points within distance epsilon from the
%                    pairs of line-segments of the dictionary atoms at each
%                    iteration of the algorithm
%
%   TODO:
%   1) Modify this function to return the sparse code learned in addition to
%   the dictionary.
%
% P = U*X
% P = d*n matrix, U = d*k matrix, X = k*n matrix
%

[d, n] = size(P);

if nargin < 3
    stopping_func = @max;
end

if nargin < 2
    % Choosing epsilon
    closest = pdist2(P', P', 'euclidean', 'Smallest', 2);
    %epsilon = min(closest(2, :)); % Dist between closest two points
    epsilon = mean(closest(2, :)); % Avg distance between pairs of closest points
end

% Learn the dictionary using the greedy distance to line segments algorithm
U = zeros(d, n);
dist_array = zeros(1, n);
avg_dist_array = zeros(1, n);
count_inactive = zeros(1, n);

r = randi([1, n], 1, 1);

U(:, 1) = P(:, r);
P(:, r) = [];
D = compute_dist_closest_line(U(:, 1), P); % Set of distances from each point in P to the closest line-seg in U
[max_dist, max_index] = max(D);
dist_array(1) = max_dist;
avg_dist_array(1) = mean(D);
count_inactive(1) = sum(D <= epsilon) + 1;

flag = 0;
for i = 2:n
    timer = tic;
    %if max_dist <= epsilon
    if stopping_func(D) <= epsilon
        flag = 1;
        break
    end
    
    U(:, i) = P(:, max_index);
    P(:, max_index) = [];
    
    if i == n
        dist_array(i) = 0;
        avg_dist_array(i) = 0;
        count_inactive(i) = n;
    end
    % Update the distance vector D due to the newly generated line segments
    % in U
    D(max_index) = [];
    dist = Inf*ones(1, n-i);
    for j = (i-1):-1:1
        [~, temp_dist] = compute_dist_point_to_line(P, U(:, j), U(:, i));     
        dist = min(dist, temp_dist);
    end
    D = min(D, dist);
    [max_dist, max_index] = max(D);
    dist_array(i) = max_dist;
    avg_dist_array(i) = mean(D);
    count_inactive(i) = sum(D <= epsilon) + i;
    
    fprintf('End of iteration:%d\n', i);
    toc(timer);
end

% for i = 2:n
%     tic;
%     if max_dist <= epsilon
%         flag = 1;
%         break
%     end
%     
%     U(:, i) = P(:, max_index);
%     P(:, max_index) = [];
%     
%     if i == n
%         dist_array(i) = 0;
%         avg_dist_array(i) = 0;
%         count_inactive(i) = n;
%     end
%     
%     D = compute_dist_closest_line(U(:, 1:i), P);
%     [max_dist, max_index] = max(D);
%     dist_array(i) = max_dist;
%     avg_dist_array(i) = mean(D);
%     count_inactive(i) = sum(D <= epsilon) + i;
%     
%     fprintf('End of iteration:%d\n', i);
%     toc;
% end

if flag == 1
    U = U(:, 1:(i-1));
    dist_array  = dist_array(1:(i-1));
    avg_dist_array  = avg_dist_array(1:(i-1));
    count_inactive = count_inactive(1:(i-1));
end
end


