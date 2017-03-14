function [recon] = visualize_reconstruction(D, Y, algorithm_id, iterations)
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
%   iterations   - #iterations of the sparse coding algorithm to be run
%   
%   OUTPUT:
%   recon - the reconstruction of each point from Y for each iteration of
%           the sparse coding algorithm (size = size(Y)*#iterations)
%




end

