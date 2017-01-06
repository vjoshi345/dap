% Compute the singular value decomposition of the input matrix and then
% reduce the number of singular vectors iteratively. At each stage we plot
% the forbenius norm of the difference between the original matrix and its
% reconstruction

% P - input matrix (mxn), m < n
% Full SVD - P = U*S*V' where U is (mxm), S is (mxn) and V is (nxn)
clear all;

% Importing data and converting to the matrix form
P = csvread('ionosphere_mod.csv');
P = P'; 
[d, n] = size(P);

% Full SVD
[U, S, V] = svd(P);
C_full = norm((P - U*S*V'), 'fro');
fprintf('Norm for full SVD: %3.3f\n', C_full);

% Economy size SVD
[U, S, V] = svd(P, 'econ');
C_econ = norm((P - U*S*V'), 'fro');
fprintf('Norm for economy size SVD: %3.3f\n', C_econ);

% Reduced SVD
cost = zeros(1, d);
cost(d) = C_econ;

for i = 1:(d-1)
    U(:, (d - i + 1)) = [];
    V(:, (d - i + 1)) = [];
    S(:, (d - i + 1)) = [];
    S((d - i + 1), :) = [];
    
    cost(d-i) = norm((P - U*S*V'), 'fro');
end

plot(cost);
title(['Reconstruction error (Frobenius norm), d = ' num2str(d) ', n = ' num2str(n)]);
xlabel('No. of singular vectors');
ylabel('Reconstruction error');
    
    
    
    
    
    
    
    
    
    
    
    
    