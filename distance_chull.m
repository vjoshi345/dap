function dist = distance_chull(P, q, m)
%distance_chull - Calculates the distance of a set of query points to the
%convex hull of a set of points
% P = matrix of the set of points (column vectors) forming the convex hull
% q = the query point (column vector) or query matrix (set of column
% vectors)
% m = no. of iterations to run
% dist = the vector of evaluated distances 

% Initial condition - find points t in P which are closest to query points q
%[d, n] = size(P); % Number of [dimensions, points]
[~, n] = size(q); % No. of query points

[~, min_index] = pdist2(P', q', 'euclidean', 'Smallest', 1);
t = P(:, min_index); % Matrix of closest points

% for i = 1:n
%     temp_dist = norm(q - P(:, i));
%     if temp_dist < min_dist
%         min_dist = temp_dist;
%         t_0 = P(:, i);
%     end
% end

% Compute the closest point from q to the convex hull of P
dist = zeros(1, n); % Vector of distances from q to conv(P)
for i = 1:m
    v = q - t;
    %v = v/norm(v); 
    v = normc(v);
    
    %max_value = -Inf;
    %p = zeros(d, 1);
%     for j = 1:n
%         temp_value = dot(v, P(:, j));
%         if temp_value > max_value
%             max_value = temp_value;
%             p = P(:, j);
%         end
%     end
    [~, max_index] = max(v'*P, [], 2); % This is a column vector
    p = P(:, max_index');
    for j = 1:n
        [t(:, j), dist(j)] = point_to_line(q(:, j), t(:, j), p(:, j));
    end
    %[t, dist] = point_to_line(q, t, p);
end
end

