function [error] = dl_nn_classif(data_path, labels_path, dl_algo, param)
% DL_NN_CLASSIF Dictionary learning based nearest neighbour classification.
%   
%   INPUT:
%   data    - classpath to the input dataset as a matrix with columns as
%             datapoints and rows as features
%   labels  - classpath to the class labels for the points
%   dl_algo - algorithm used for dictionary learning
%             1-dp, 2-dl, 3-dch, 4-dchperceptron
%   param   - number of optional arguments as a struct
%       stopping_func - (optional) either @max or @mean - to be used as the
%                       stopping criterion (default = @max)
%       nneigbours    - (optional)#neigbours to take majority vote on 
%                       (default=1). If nneighbours = -1, then
%                       takes majority vote on the 10 atoms in the sparse
%                       coding of the datapoint
%   
%   OUTPUT:
%   error - classification error on the data
%

rng('default');

Y = csvread(data_path);
[d, n] = size(Y);
labels = csvread(labels_path);
assert(size(labels, 1) == n, '#labels does not equal #points!');

[~, file_name, ~] = fileparts(data_path);
fprintf('\n');
fprintf('-------------------------------------------------------------\n');
display(['Dataset name:' file_name]);
display(['No. of points(n):' num2str(n)]);
display(['No. of dimensions(d):' num2str(d)]);

closest = pdist2(Y', Y', 'euclidean', 'Smallest', 2);
epsilon = mean(closest(2, :)); % Avg distance between pairs of closest points
display(['Error tolerance (epsilon)=' num2str(epsilon)]);

if (~isfield(param, 'stopping_func'))
    stopping_func = @max;
else
    stopping_func = param.stopping_func;
end
display(['Stopping function=' func2str(stopping_func)]);

if (~isfield(param, 'nneighbours'))
    nneighbours = 1;
else    
    nneighbours = param.nneighbours;
    if nneighbours == -1
        assert(dl_algo == 3 || dl_algo == 4, 'For variable nneighbours need to use either dch or dchperceptron!');        
        if dl_algo == 3
            method = @compute_dist_chull;
        else
            method = @compute_dist_chull_perceptron;
        end
    end
end
display(['# neighbours=' num2str(nneighbours)]);

%--------- DICTIONARY LEARNING ---------------
algorithm_list = {@dp, @dl, @dch, @dch};
algorithm = algorithm_list{dl_algo};
switch dl_algo
    case 1
        algorithm_name = [func2str(algorithm) '-' func2str(stopping_func)];
        display(['Algorithm chosen:' algorithm_name]);
        [selected, ~, ~, ~, ~] = algorithm(Y, epsilon, stopping_func);        
        D = Y(:, selected);
    case 2
        algorithm_name = [func2str(algorithm) '-' func2str(stopping_func)];
        display(['Algorithm chosen:' algorithm_name]);
        [D, ~, ~, ~] = algorithm(Y, epsilon, stopping_func);        
    case 3
        algorithm_name = [func2str(algorithm) '-' func2str(stopping_func)];
        display(['Algorithm chosen:' algorithm_name]);
        [D, ~, ~, ~] = algorithm(Y, 1, epsilon, 10, stopping_func);        
    case 4
        algorithm_name = [func2str(algorithm) 'perceptron-' func2str(stopping_func)];
        display(['Algorithm chosen:' algorithm_name]);
        [D, ~, ~, ~] = algorithm(Y, 2, epsilon, 10, stopping_func);
    otherwise
        disp('Incorrect input');
        return
end
k = size(D, 2);
disp(['Size of dictionary learned:' num2str(k)]);
% 
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist2-deskewed_Dictionary_dl-mean.mat', 'D');

[~, dict_idx] = ismember(D', Y', 'rows');
unselected = setdiff(1:n, dict_idx);
X = Y(:, unselected);
labels_D = labels(dict_idx);
labels_X = labels(unselected);

%----------- NN-CLASSIFICATION ------------------
if nneighbours >= 1
    [~, pred_idx] = pdist2(D', X', 'euclidean', 'Smallest', nneighbours);
    pred = zeros(n-k, 1);
    for i = 1:(n-k)
        index = pred_idx(:, i);
        index = index(index > 0);
        pred(i) = mode(labels_D(index));
    end
    error = sum(~(pred == labels_X))*100/(n-k);
else
    count_atoms = ones(n-k, 1);
    recon = zeros(d, n-k);
    covered_id = [];
    for j = 1:10            
        [dist, t, atom_idx] = method(D, X, j);
        if isequal(stopping_func, @max)
            threshold = epsilon;
        else
            threshold = 2*epsilon - max(dist);
        end
        threshold = max(threshold, 0);
        if j < 10
            count_atoms(dist > threshold) = count_atoms(dist > threshold) + 1;
        end
        new_id = setdiff(find(dist <= threshold), covered_id);
        recon(:, new_id) = t(:, new_id);
        covered_id = [covered_id new_id];
    end
    
    %[~, ~, atom_idx] = method(D, X, 10);
%     pred = zeros(n-k, 1);
%     for i = 1:(n-k)
%         index = atom_idx(:, i);
%         index = index(index > 0);
%         max_id = min(count_atoms(i), size(index, 1));
%         index = index(1:max_id);
%         pred(i) = mode(labels_D(index));
%     end
    [~, pred_idx] = pdist2(D', recon', 'euclidean', 'Smallest', 1);
    pred = labels_D(pred_idx);
    error = sum(~(pred == labels_X))*100/(n-k);
end

results = [n, d, epsilon, k, nneighbours, error];
string = sprintf('%0.5f,', results);
string = [file_name ',' algorithm_name ',' string '\n'];
fid = fopen('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\results12.csv', 'a');
fprintf(fid, string);
fclose(fid);    














end

