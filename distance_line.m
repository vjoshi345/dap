function dist = distance_line(P, q)
%distance_line - Calculates the farthest distance from a given query point  
%to the set of line segments formed by a set of points (for n points "n 
%choose 2" set of line segments)  
% P = matrix of the set of points (column vectors) 
% q = the query point (column vector)
% dist = the evaluated distance

max_dist = 0;
[d, n] = size(P); % Number of [dimensions, points]
t_0 = zeros(d, 1); % Farthest point

for i = 1:(n-1)
	for j = (i+1):n
		[temp_point, temp_dist] = point_to_line(q, P(:, i), P(:, j));
		if temp_dist > max_dist
			max_dist = temp_dist
			t_0 = temp_point
		end
    end
end

dist = max_dist;
end

