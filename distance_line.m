function [dist, index] = distance_line(U, P)
%distance_line - Given matrices U and P, this function gives the largest 
%distance from points in P to the closest line segment induced by pairs of
%points in U.  
% U = matrix of the set of points (column vectors)
% P = matrix of the set of points (column vectors) 
% dist = the evaluated distance
% index = index of the point in P which is farthest from the closest line 
% segments in U
tic;
[~, k] = size(U); % Number of [dimensions, points]
n = size(P, 2); % Number of points

closest_lineseg_dist = zeros(1, n);

for a = 1:n
    min_dist = Inf;
    for b = 1:(k-1)
        for c = (b+1):k
            [~, temp_dist] = point_to_line(P(:, a), U(:, b), U(:, c));
            if temp_dist < min_dist
                min_dist = temp_dist;
            end
        end
    end
    closest_lineseg_dist(a) = min_dist;
end

[dist, index] = max(closest_lineseg_dist);
toc;
end


