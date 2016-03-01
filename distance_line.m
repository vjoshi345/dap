function [dist, index] = distance_line(U, P)
%distance_line - Given matrices U and P, this function gives the largest 
%distance from points in P to the closest line segment induced by pairs of
%points in U.  
% U = matrix of the set of points (column vectors)
% P = matrix of the set of points (column vectors) 
% dist = the evaluated distance
% index = index of the point in P which is farthest from the closest line 
% segments in U

[~, k] = size(U); % Number of [dimensions, points]
n = size(P, 2); % Number of points

if k == 1
    [~, temp] = point_to_line(P, U, U);
    [dist, index] = max(temp);
else
    lin_seg_count = k*(k-1)/2;
    lineseg_dist = zeros(lin_seg_count, n); % Distance from each point to 
    i = 0;                                  % each line segment
    for b = 1:(k-1)
        for c = (b+1):k
            i = i + 1;
            [~, lineseg_dist(i, :)] = point_to_line(P, U(:, b), U(:, c));
        end
    end

    min_lineseg_dist = min(lineseg_dist, [], 1); % Vector of closest 
                                                 % distances for each point
    [dist, index] = max(min_lineseg_dist);
end
end


