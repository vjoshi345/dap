function [recon] = visualize_reconstruction(D, Y, algorithm_id, param)
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
%   param        - structure that includes a variety of optional arguments
%       niter    - (optional) max #iterations of the sparse coding 
%                  algorithm to be run
%       epsilon  - (optional) error threshold for the reconstruction
%   
%   OUTPUT:
%   recon - the reconstruction of each point from Y for each iteration of
%           the sparse coding algorithm (size = size(Y)*#iterations)
%

% Check that the input arguments are correctly specified
[d, n] = size(Y);
assert(d == size(D, 1), 'Dimensions of dictionary and original points do not match!');

assert(algorithm_id == 1 || algorithm_id == 2, 'algorithm_id should be either 1 or 2!');

if(~isfield(param, 'niter'))
    niter = 10;
else
    niter = param.niter;
end

if(~isfield(param, 'epsilon'))
    epsilon = 0;
else
    epsilon = param.epsilon;
end

method_list = {@compute_dist_chull, @compute_dist_chull_perceptron};
method = method_list{algorithm_id};
algorithm_list = {'dch', 'dchperceptron'};
algorithm_name = algorithm_list{algorithm_id};
shape = sqrt(d);

% Get the reconstructed points for each level of sparsity
recon = zeros(d, n, niter);
dist = zeros(niter, n);
idx = false(niter, n);
idx(1, :) = true;
for i = 1:niter
     [dist(i, :), recon(:, :, i)] = method(D, Y, i);
     idx(i, :) = dist(i, :) > epsilon;
end

% Plot the reconstructed points alongwith the original
for i = 1:n
    recon_fig = figure('visible', 'off');
    title(['Reconstruction of a point using ' algorithm_name]);
    curriter = sum(idx(:, i));
    if curriter == 0
        curriter = 1;
    end
    subplot(2, curriter, 1:curriter);
    curr = reshape(Y(:, i), [shape, shape]);
    imshow(curr);
    %imshow(ifft2(curr));
    title('Original point');
    for j = 1:curriter
        subplot(2, curriter, curriter+j);
        curr = reshape(recon(:, i, j), [shape, shape]);
        imshow(curr);
        %imshow(ifft2(curr));
        title(['Sp:' num2str(j)]);
    end
    saveas(recon_fig, ['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\to_publish\Visualizations\Reconstruction\Deskewed\mnist1-deskewed_dch_reconstruction' num2str(i) '.jpg']);
end




























