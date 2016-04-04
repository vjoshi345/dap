function dist = distance_chull_perceptron(P, q, m)
%distance_chull_perceptron - Calculates the distance of the set of query 
%points q to the convex hull of the set of points P using the perceptron
%method
% P = matrix of the set of points (column vectors) forming the convex hull
% q = the query point (column vector) or query matrix (set of column
% vectors)
% m = no. of iterations to run
% dist = the vector of evaluated distances

% Initial condition - find points t in P which are closest to query points q
[~, n] = size(q); % No. of query points
[~, min_index] = pdist2(P', q', 'euclidean', 'Smallest', 1);
t = P(:, min_index); % Matrix of closest points

%dist = zeros(1, n);

for i = 1:(m-1)
    temp = i*t;
    diff = temp - i*q;
    [~, min_index] = pdist2(-P', diff', 'euclidean', 'Smallest', 1);
    
    points = P(:, min_index);
    t = (temp + points)/(i+1);
end

dist = sqrt(sum((t-q).^2, 1));
end



