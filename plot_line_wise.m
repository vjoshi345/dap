% P = U*X
% P = d*n matrix, U = d*k matrix, X = k*n matrix
% This program plots the max distance from the subset of points (U) chosen
% using three different algorithms as a function of the size of the subset.

tic;
% Setting Constants
epsilon = 1;

% Importing data and converting to the matrix form
P = csvread('ionosphere_mod.csv');
P = P'; 
[d, n] = size(P);

Q = P; % We will manipuate P and keep a copy of it in Q for later use

%%%%%%%%%%%%%%%%%% ALGORITHM-2 Line-wise max distance %%%%%%%%%%%%%%%%%%%%
% Computing columns of U 
r = randi([1, n], 1, 1);
U = zeros(d, n);
U(:, 1) = P(:, r);
P = P(:, [(1:r-1), (r+1:end)]);

dist_array = zeros(1, n);
dist_array(1) = distance_point(P, U(:, 1));

i = 1;
while i < n
    [D, I] = pdist2((U(:, 1:i))', P', 'euclidean', 'Smallest', 1);
    [max_dist, max_index] = max(D); % Max distance among the closest points
                                    % of U from P
     
    if max_dist <= epsilon
        break
    else
        U(:, i+1) = P(:, max_index);
        P = P(:, [(1:max_index-1), (max_index+1:end)]);
        dist_array(i+1) = max_dist;
    end
    fprintf('End of iteration:%d\n', i);
    i = i + 1;
end
U = U(:, 1:i);
dist_array = dist_array(1:i);

%%% Plotting the distance as a function of the size of subset chosen
plot(dist_array);

% %%%%%%%%%%%%%%%%%% ALGORITHM-2 Line segment-wise max distance %%%%%%%%%%%%%%%%%%%%
% P = Q; % We will manipuate P and keep a copy of it in Q for later use
% % Computing columns of U 
% r = randi([1, n], 1, 1);
% U = zeros(d, n);
% U(:, 1) = P(:, r);
% P = P(:, [(1:r-1) (r+1:end)]);
% 
% line_dist_array = zeros(1, n);
% line_dist_array(1) = distance_line(P, U(:, 1));
% 
% max_dist = -Inf;
% far_point = zeros(d, 1);
% i = 0;
% while true
%     i = i + 1;
%     max_dist = distance_line(U(:, 1:i), P(:, 1));
%     far_point = P(:, 1);
% 	n = size(P, 2);
%     for j = 2:n
%         temp_dist = distance_line(U(:, 1:i), P(:, j));
%         if temp_dist > max_dist
%             max_dist = temp_dist;
%             far_point = P(:, j);
% 			f = j;
%         end
%     end
%     if max_dist <= epsilon
%         break
%     else
%         U(:, i+1) = far_point;
% 		line_dist_array(i+1) = max_dist;
% 		P = P(:, [(1:f-1) (f+1:end)]);
%     end
%     fprintf('End of iteration:%d\n', i);
% end
% U = U(:, 1:i);
% line_dist_array = line_dist_array(1:i);

% %%% Plotting the distance as a function of the size of subset chosen
% plot(line_dist_array);
    
% % Evaluate X using matrix multiplication
% c = 0.1;
% U = U + c*eye(size(U));
% %X = (U'*U)\(U'*P);
% X = U\P;

% % Cost (Frobenius norm) ||P - U*X||^2 
% C = norm((P - U*X), 'fro');

% % Sparsity
% nonzero_count = sum(sum(abs(X) > 1e-5));
% sparsity_coeff = nonzero_count/numel(X);

% % Output
% fprintf('No. of points(n): %d\t', n);
% fprintf('No. of dimensions(d): %d\n', d);
% fprintf('No. of clusters(k): %d\n', size(U, 2));
% fprintf('Cost ||P - UX||: %3.3f\n', C);
% fprintf('Sparsity: %3.3f\n', sparsity_coeff);

toc;