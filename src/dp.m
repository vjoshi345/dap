function [selected, sparse_code, dist_array, avg_dist_array, count_inactive] = dp(P, epsilon, stopping_func)
% DP Learns a dictionary for an input dataset using the distance from
% points algorithm.
%
%   INPUT:
%   P             - the input dataset as a matrix with columns as 
%                   datapoints and rows as dimensions
%   epsilon       - error tolerance for each datapoint
%   stopping_func - either @max or @mean which is to be used as the
%                   stopping criterion (default = @max)
%
%   OUTPUT:
%   selected       - the index of atoms selected in the dictionary (vector)
%   sparse_code    - vector specifying sparse code for the data points
%   dist_array     - distance of farthest point from dictionary atoms at
%                    each iteration of the algorithm
%   avg_dist_array - average distance of points from dictionary atoms at
%                    each iteration of the algorithm
%   count_inactive - number of points within distance epsilon from
%                    dictionary atoms at each iteration of the algorithm
%
% P = U*X
% P = d*n matrix, U = d*k matrix, X = k*n matrix
%

[~, n] = size(P);

if nargin < 3
    stopping_func = @max;
end

% Choose epsilon if argument is missing 
if nargin < 2
    closest = pdist2(P', P', 'euclidean', 'Smallest', 2);
    %epsilon = min(closest(2, :)); % Dist between closest two points
    epsilon = mean(closest(2, :)); % Avg dist between closest pair of points
end

% Learn the dictionary and sparse code simultaneously
selected = false(1, n);
not_selected = true(1, n);
not_selected_index_array = 1:n;
sparse_code = 1:n;

r = randi([1, n], 1, 1);
selected(r) = true;
not_selected(r) = false;
not_selected_index_array(r) = [];

dist_array = zeros(1, n);
avg_dist_array = zeros(1, n);
count_inactive = zeros(1, n);

% I contains indices from selected. max_index is one of the indices from
% not_selected.
[D, I] = pdist2(P(:, selected)', P(:, not_selected)', 'euclidean', 'Smallest', 1);
[max_dist, max_index] = max(D);

dist_array(1) = max_dist;
avg_dist_array(1) = mean(D);
count_inactive(1) = sum(D <= epsilon) + 1;

flag = 0;
for i = 2:n
    %if max_dist <= epsilon
    if stopping_func(D) <= epsilon
        flag = 1;
        break
    end
    
    selected(not_selected_index_array(max_index)) = true;
    not_selected(not_selected_index_array(max_index)) = false;
    %not_selected_index_array(max_index) = [];
    
    if i == n
        dist_array(i) = 0;
        avg_dist_array(i) = 0;
        count_inactive(i) = n;
        break
    end
    
    % Update the distance vector D and index vector I due to the newly 
    % selected point in the dictionary
    %[D, I] = pdist2(P(:, selected)', P(:, not_selected)', 'euclidean', 'Smallest', 1);
    D(max_index) = [];
    I(max_index) = [];
    [dist, ~] = pdist2(P(:, not_selected_index_array(max_index))', P(:, not_selected)', 'euclidean', 'Smallest', 1);
    not_selected_index_array(max_index) = [];
    for j = 1:size(dist, 2)
        if dist(j) < D(j)
            I(j) = i;
            D(j) = dist(j);
        end
    end    
    [max_dist, max_index] = max(D);
    dist_array(i) = max_dist;
    avg_dist_array(i) = mean(D);
    count_inactive(i) = sum(D <= epsilon) + i;
    
    fprintf('End of iteration:%d\n', i);
end

if flag == 1
    dist_array  = dist_array(1:(i-1));
    avg_dist_array  = avg_dist_array(1:(i-1));
    count_inactive = count_inactive(1:(i-1));
end

% Generate the sparse code
selected_index_array = 1:n;
selected_index_array = selected_index_array(selected);
sparse_code(not_selected) = selected_index_array(I);

end



