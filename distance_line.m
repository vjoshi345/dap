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

lin_seg_count = k*(k-1)/2;
lineseg_dist = zeros(lin_seg_count, n);

i = 0;
for b = 1:(k-1)
    for c = (b+1):k
        i = i + 1;
        [~, lineseg_dist(i, :)] = point_to_line(P, U(:, b), U(:, c));
    end
end

min_lineseg_dist = min(lineseg_dist);

[dist, index] = max(min_lineseg_dist);
toc;
end


