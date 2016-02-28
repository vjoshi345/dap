% P = U*X
% P = d*n matrix, U = d*k matrix, X = k*n matrix
tic;
% Setting Constants
epsilon = 1;

% Importing data and converting to the matrix form
P = csvread('ionosphere_mod.csv');
P = P'; 
%P = P(:, 1:50);
[d, n] = size(P);

% Computing columns of U (greedy algorithm)
r = randi([1, n], 1, 1);
U = zeros(d, n);
U(:, 1) = P(:, r);
P(:, r) = [];

dist_array = zeros(1, n);
dist_array(1) = distance_chull;
max_dist = -Inf;
far_point = zeros(d, 1);
i = 0;


while true
    i = i + 1;
    max_dist = distance_chull(U, P(:, 1), 100);
    dist_array(1) = max_dist;
    far_point = P(:, 1);
    for j = 2:n
        temp_dist = distance_chull(U, P(:, j), 100);
        if temp_dist > max_dist
            max_dist = temp_dist;
            far_point = P(:, j);
        end
    end
    if max_dist <= epsilon
        break
    else
        U(:, i+1) = far_point;
        dist_array(i+1) = max_dist;
    end
    fprintf('End of iteration:%d\n', i);
end
U = U(:, 1:i);
dist_array = dist_array(1:i);

plot(dist_array);
    
% % Evaluate X using matrix multiplication
% c = 0.1;
% U = U + c*eye(size(U));
% %X = (U'*U)\(U'*P);
% X = U\P;
% 
% % Cost (Frobenius norm) ||P - U*X||^2 
% C = norm((P - U*X), 'fro');
% 
% % Sparsity
% nonzero_count = sum(sum(abs(X) > 1e-5));
% sparsity_coeff = nonzero_count/numel(X);
% 
% % Output
% fprintf('No. of points(n): %d\t', n);
% fprintf('No. of dimensions(d): %d\n', d);
% fprintf('No. of clusters(k): %d\n', size(U, 2));
% fprintf('Cost ||P - UX||: %3.3f\n', C);
% fprintf('Sparsity: %3.3f\n', sparsity_coeff);

toc;