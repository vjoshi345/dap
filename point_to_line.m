function [p, d] = point_to_line(q, a, b)
%point_to_line Returns the point (and its distance) which is closest to a
%query point on a line segment
%   q = query point
%   a = end point of line segment
%   b = end point of line segment
%   p = point closest to q on the line segment ab
%   d = distance from q to p

if norm(b-a) == 0
    p = a;
    d = norm(q-a);
else
    x = q - a;
    y = b - a;

    t = dot(x, y)/dot(y, y);

    p = a + t*y;

    d = norm(x - t*y);
end
end

