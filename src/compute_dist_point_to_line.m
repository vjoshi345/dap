function [p, dist] = compute_dist_point_to_line(q, a, b)
% COMPUTE_DIST_POINT_TO_LINE Calculates the distance from a query point (or
% a set of points) to a line segment (or a set of line segments). Also 
% returns the closest point to query on the line segment.
%
%    INUPT:
%    q - query point (or matrix of points) (dxn)
%    a - end point of line segment (matrix dxn)
%    b - end point of line segment (matrix dxn)
%    
%    OUTPUT:
%    p    - point on ab which is closest to q (matrix dxn)
%    dist - distance from q to to line segment ab (vector 1xn)
%
%   TODO:
%   Modify this function to work with q, a, and b all matrices. Each column
%   of q is a query point. Corresponding columns of a and b are end points
%   of line segements corresponding to a point in q.

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
