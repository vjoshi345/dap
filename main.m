% P = U*X
% P = d*n matrix, U = d*k matrix, X = k*n matrix

% Setting Constants
epsilon = 0.01;

% Importing data and converting to the matrix form
P = csvread('ionosphere_mod.csv');
P = P';

% Computing columns of U (greedy algorithm)

% Evaluate X using matrix multiplication


% Cost (Frobenius norm) ||P - U*X||^2

% Output
