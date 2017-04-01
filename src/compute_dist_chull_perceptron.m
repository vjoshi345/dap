function [dist, t, atom_idx] = compute_dist_chull_perceptron(P, q, niter)
% COMPUTE_DIST_CHULL_PERCEPTRON - Calculates the distance of a set of query 
% points to the convex hull of a set of points using the perceptron method.
%
%    INPUT: 
%    P     - matrix of the set of points (column vectors) forming the 
%            convex hull
%    q     - matrix of query points represented as column vectors
%    niter - no. of iterations to get approximate distance to convex hull
%            N.B.: note that niter here reflects the sparsity of
%            reconstruction
%
%    OUTPUT:
%    dist     - the vector of evaluated distances (1xn)
%    t        - the matrix of reconstructed points (dxn)
%    atom_idx - matrix of dictionary atoms (indicies) chosen during the
%               sparse coding stage (niter x #query points)
%
%   TODO:
%   Need to check different possible interpretations of the perceptron
%   algorithm and see which one to use here (right now we are using one
%   of the variations).
%

[~, n] = size(q);
atom_idx = zeros(niter, n);

% Initial condition - find points t in P which are closest to query points 
[~, min_index] = pdist2(P', q', 'euclidean', 'Smallest', 1);
t = P(:, min_index); % Matrix of closest points
atom_idx(1, :) = min_index; 
for i = 1:(niter-1)
    % Find the next best point to add in t
    temp = i*t;
    diff = temp - (i+1)*q;
    [~, min_index] = pdist2(-P', diff', 'euclidean', 'Smallest', 1);
    points = P(:, min_index);
    
    % Checking whether adding "points" leads to reduced distance to q
    curr_dist = sqrt(sum((t-q).^2, 1));
    new_t = (temp + points)/(i+1);
    new_dist = sqrt(sum((new_t-q).^2, 1));
    qual_index = new_dist < curr_dist;
    
    % Adding the qualified points in t 
    t(:, qual_index) = (temp(:, qual_index) + points(:, qual_index))/(i + 1);
    atom_idx(i+1, :) = qual_index;
end
dist = sqrt(sum((t-q).^2, 1));
end



