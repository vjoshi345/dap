function [recon] = visualize_reconstruction(D, Y, algorithm_id, niter)
% VISUALIZE_RECONSTRUCTION Plots the reconstruction of a point at each
% iteration using a fixed dictionary and a sparse coding algorithm.
% N.B.: this function is meaningful only for an image dataset.
%
%   INPUT:
%   D            - the input dictionary(rows = features, columns = atoms)
%   Y            - data matrix which will be reconstructed(columns =
%                  points, rows = features)
%   algorithm_id - algorithm used for sparse coding
%                  1-dch, 2-dchperceptron
%   niter        - #iterations of the sparse coding algorithm to be run
%   
%   OUTPUT:
%   recon - the reconstruction of each point from Y for each iteration of
%           the sparse coding algorithm (size = size(Y)*#iterations)
%

assert(algorithm_id == 1 || algorithm_id == 2, 'algorithm_id should be either 1 or 2!');
[d, n] = size(Y);
assert(d == size(D, 1), 'Dimensions of dictionary and original points do not match!');

method_list = {@compute_dist_chull, @compute_dist_chull_perceptron};
method = method_list{algorithm_id};
algorithm_list = {'dch', 'dchperceptron'};
algorithm_name = algorithm_list{algorithm_id};

% Get the reconstructed points for each level of sparsity
recon = zeros(d, n, niter);
for i = 1:niter
     [~, recon(:, :, i)] = method(D, Y, i);
end

% Plot the reconstructed points alongwith the original
for i = 1:n
    figure();
    subplot(2, niter, 1:niter);
    curr = reshape(Y(:, i), [28, 28]);
    imshow(curr);
    title('Original point');
    for j = 1:niter
        subplot(2, niter, niter+j);
        curr = reshape(recon(:, i, j), [28, 28]);
        imshow(curr);
        title(['Sparsity:' num2str(j)]);
    end
    suptitle(['Reconstruction of a point using ' algorithm_name]);
end




























