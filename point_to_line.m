function [dist] = point_to_line(q, a, b)
%point_to_line Returns the distance from a query point (or set of points) 
%to a line segment
%   q = query point (or matrix of points) (dxn)
%   a = end point of line segment (column vector dx1)
%   b = end point of line segment (column vector dx1)
%   d = distance from q to to line segment ab

if isequal(a, b)
    %p = a;
    dist = pdist2(a', q'); % dist is (1xn)
else
    x = bsxfun(@minus, q, a); % x is (dxn)
    y = b - a; % y is (dx1)

    t = (x'*y)/(y'*y); % t is (nx1)

    %p = bsxfun(@plus, (y*t'), a); % p is (dxn)

    dist = sqrt(sum((x - (y*t')).^2, 1)); % dist is (1xn)
end
end

