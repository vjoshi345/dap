function [error] = random_nn_classif(data_path, labels_path, k, nneighbours)
% RANDOM_NN_CLASSIF Random dictionary based nearest neighbour classification.
%   
%   INPUT:
%   data_path   - classpath to the input dataset as a matrix with columns as
%                  datapoints and rows as features
%   labels_path - classpath to the class labels for the points
%   k           - the random subset size to be taken as the dictionary
%   nneighbours - #neighbours to take majority vote on
%   
%   OUTPUT:
%   error - classification error on the data
%
rng('default');

Y = csvread(data_path);
[d, n] = size(Y);
labels = csvread(labels_path);
assert(size(labels, 1) == n, '#labels does not equal #points!');
assert(k <= n, 'Dictionary size cannot be greater than #datapoints!');

assert(nneighbours >= 1, '# neighbours for kNN should be >=1');

[~, file_name, ~] = fileparts(data_path);
fprintf('\n');
fprintf('-------------------------------------------------------------\n');
display(['Dataset name:' file_name]);
display(['No. of points(n):' num2str(n)]);
display(['No. of dimensions(d):' num2str(d)]);
display(['# neighbours=' num2str(nneighbours)]);
disp(['Size of dictionary:' num2str(k)]);

%--------- DICTIONARY LEARNING ---------------
D = datasample(Y, k, 2, 'Replace', false);
[~, dict_idx] = ismember(D', Y', 'rows');
unselected = setdiff(1:n, dict_idx);
X = Y(:, unselected);
labels_D = labels(dict_idx);
labels_X = labels(unselected);

%----------- NN-CLASSIFICATION ------------------
[~, pred_idx] = pdist2(D', X', 'euclidean', 'Smallest', nneighbours);
pred = zeros(n-k, 1);
for i = 1:(n-k)
    index = pred_idx(:, i);
    % Weighting by distance
    weights = sqrt(sum(bsxfun(@minus, X(:, i), D(:, index)).^2, 1))';
    weights = 1./weights;
    weights = weights/sum(weights);
    [pot_labels, ~, ic] = unique(labels_D(index));
    weights_by_labels = accumarray(ic, weights, size(pot_labels), @(x) sum(x));
    [~, idx_label] = max(weights_by_labels);
    pred(i) = pot_labels(idx_label); 
%     pred(i) = mode(labels_D(index));
end
error = sum(~(pred == labels_X))*100/(n-k);

%------------- STORE RESULTS -----------------------
epsilon = 4.4085;
results = [n, d, epsilon, k, nneighbours, error];
string = sprintf('%0.5f,', results);
string = [file_name ',' 'randomized,' string '\n'];
fid = fopen('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\results15.csv', 'a');
fprintf(fid, string);
fclose(fid);    

end

