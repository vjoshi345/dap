function [p, dist] = compute_dist_point_to_line(q, a, b)
% COMPUTE_DIST_POINT_TO_LINE Calculates the distance from a query point (or
% a set of points) to a line segment (or a set of line segments). Also 
% returns the closest point to query on the line segment.
%
%    There are two ways to use this function:
%    1) q is a matrix and a and b are vectors. This function returns the
%    set of distances for each point in q to the line segment ab.
%    2) q is a matrix and a and b are matrices (q, a, and b have the same
%    size). This function returns the distance of each point in q to the
%    corresponding line segment from ab. For example, dist(i) = distance
%    from q(:, i) to line segment a(:, i)b(:, i).
%
%    INUPT:
%    q - query point (or matrix of points) (dxn)
%    a - end point of line segment (vector dx1)
%    b - end point of line segment (vector dx1)
%    
%    OUTPUT:
%    p    - point on ab which is closest to q (matrix dxn)
%    dist - distance from q to to line segment ab (vector 1xn)
%

if ~isvector(a) && ~isvector(b)
    x = q - a; % x is (dxn)
    y = b - a; % y is (dxn)
    t = sum(x.*y); % t is (1xn)
    den = sum(y.^2); % den is (1xn)
    nonzero_index = find(den);
    t(nonzero_index) = t(nonzero_index)./den(nonzero_index);
    p = a + bsxfun(@times, y, t); % p is (dxn)
    
    dist = sqrt(sum((q - p).^2)); % dist is (1xn)
else
    if isequal(a, b)
        p = a;
        dist = pdist2(a', q'); % dist is (1xn)
    else
        x = bsxfun(@minus, q, a); % x is (dxn)
        y = b - a; % y is (dx1)

        t = (x'*y)/(y'*y); % t is (nx1)

        p = bsxfun(@plus, (y*t'), a); % p is (dxn)

        dist = sqrt(sum((x - (y*t')).^2, 1)); % dist is (1xn)
    end
end
end

