function dist = distance_chull(P, q, m)
%distance_chull - Calculates the distance of a given query point to the
%convex hull of a set of points
% P = matrix of the set of points (column vectors) forming the convex hull
% q = the query point (column vector)
% m = no. of iterations to run
% dist = the evaluated distance

% Initial condition - find closest point t_0
min_dist = Inf;
[d, n] = size(P); % Number of [dimensions, points]
t_0 = zeros(d, 1); % Closest point

for i = 1:n
    temp_dist = norm(q - P(:, i));
    if temp_dist < min_dist
        min_dist = temp_dist;
        t_0 = P(:, i);
    end
end

% Compute the closest point from q to the convex hull of P
t = t_0;

for i = 1:m
    v = q - t;
    v = v/norm(v); 
    
    max_value = -Inf;
    p = zeros(d, 1);
    for j = 1:n
        temp_value = dot(v, P(:, j));
        if temp_value > max_value
            max_value = temp_value;
            p = P(:, j);
        end
    end
    
    [t, dist] = point_to_line(q, t, p);
end
%display(t);
end

