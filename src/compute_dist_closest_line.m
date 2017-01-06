function [dist] = compute_dist_closest_line(U, P)
% COMPUTE_DIST_CLOSEST_LINE - Given matrices U and P, this function gives 
% the distance from points in P to the closest line segment induced by 
% pairs of points in U.
%
%    INPUT:
%    U - matrix of set of points (as column vectors) in the dictionary
%    P - matrix of set of points (as column vectors) outside of dictionary
%
%    OUTPUT:
%    dist  - the vector of evaluated distances
%    index - index of the point in P which has the maximum distance from 
%            its closest line segment
%

[~, k] = size(U); % Number of [dimensions, points]
n = size(P, 2); % Number of points

if k == 1
    [~, temp] = compute_dist_point_to_line(P, U, U);
    %[dist, index] = max(temp);
    dist = temp;
else
    lin_seg_count = k*(k-1)/2;
    lineseg_dist = zeros(lin_seg_count, n); % Distance from each point to 
    i = 0;                                  % each line segment
    for b = 1:(k-1)
        for c = (b+1):k
            i = i + 1;
            [~, lineseg_dist(i, :)] = compute_dist_point_to_line(P, U(:, b), U(:, c));
        end
    end
    
    dist = min(lineseg_dist, [], 1);
    %min_lineseg_dist = min(lineseg_dist, [], 1); % Vector of closest 
                                                 % distances for each point
    %[dist, index] = max(min_lineseg_dist);
end
end


