function [dist] = compute_dist_chull(P, q, niter)
% COMPUTE_DISTANCE_CHULL Calculates the distance of a set of query points 
% to the convex hull of a set of points.
%
%   INPUT: 
%   P     - matrix of the set of points (as column vectors) which form
%                the convex hull.
%   q     - matrix of query points represented as column vectors.
%   niter - number of iterations to get approximate distance to convex
%                hull.
%
%   OUTPUT:
%   dist - the vector of evaluated distances.
%   
%   TODO:
%   Get rid of the for loop at the end. Requires changing the compute_dist_
%   point_to_line function.

% Initial condition - find points t in P which are closest to query points 
[~, n] = size(q); % No. of query points
[~, min_index] = pdist2(P', q', 'euclidean', 'Smallest', 1);
t = P(:, min_index); % Matrix of closest points

% Compute the closest point from q to the convex hull of P
dist = zeros(1, n); % Vector of distances from q to conv(P)
for i = 1:niter
    v = q - t;
    v = normc(v);
    [~, max_index] = max(v'*P, [], 2); % This is a column vector
    p = P(:, max_index');
    for j = 1:n
        [t(:, j), dist(j)] = compute_dist_point_to_line(q(:, j), t(:, j), p(:, j));
    end
end
end

