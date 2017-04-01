function [dist, t, atom_idx] = compute_dist_chull(P, q, niter)
% COMPUTE_DISTANCE_CHULL Calculates the distance of a set of query points 
% to the convex hull of a set of points.
%
%   INPUT: 
%   P     - matrix of the set of points (as column vectors) which form the
%           convex hull
%   q     - matrix of query points represented as column vectors
%   niter - no. of iterations to get approximate distance to convex hull
%           N.B.: note that niter here reflects the sparsity of
%           reconstruction
%
%   OUTPUT:
%   dist     - the vector of evaluated distances (1xn)
%   t        - the matrix of reconstructed points (dxn)
%   atom_idx - matrix of dictionary atoms (indicies) chosen during the
%              sparse coding stage (niter x #query points)
%


% Initial condition - find points t in P which are closest to query points 
[~, n] = size(q); % No. of query points
atom_idx = zeros(niter, n);
[dist, min_index] = pdist2(P', q', 'euclidean', 'Smallest', 1);
t = P(:, min_index); % Matrix of closest points
atom_idx(1, :) = min_index;

% Compute the closest point from q to the convex hull of P
%dist = zeros(1, n); % Vector of distances from q to conv(P)
for i = 1:(niter-1)
    v = q - t;
    v = normc(v);
    [~, max_index] = max(v'*P, [], 2); % This is a column vector
    p = P(:, max_index');
    [t, dist] = compute_dist_point_to_line(q, t, p);
    atom_idx(i+1, :) = max_index';
end
end

