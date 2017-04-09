function [out] = main(data_path, algorithm_id, param)
% MAIN Function to perform sparse dictionary learning on a dataset with a 
% user specified algorithm.
%
%   INPUT:
%   data_path          - full path to the dataset stored as a .csv or .mat 
%                        file
%   algorithm_id       - algorithm to perform dictionary learning
%                        (specified as a number)
%                        1 - distance from points
%                        2 - distance from line segments
%                        3 - distance from convex hull
%                        4 - distance from convex hull (perceptron method)
%                        5 - KSVD algorithm
%                        6 - randomized dictionary learning
%   param              - structure that includes a variety of optional 
%                        arguments
%       max_sparsity       - optional argument specifying the no. of 
%                            iterations to compute the distance to convex 
%                            hull for the dch algorithms (default = #dimensions of input data)
%       stopping_criterion - optional argument specifying either the max(1)
%                            or mean(2) as the stopping criteria (default = max)
%       epsilon            - optional argument specifying the error tolerance
%       dictionary_size    - (optional) if KSVD or randomized algorithm is used
%                            need to specify a dictionary size (default = 50)
%       sparse_algo        - (optional) if randomized dictionary learning
%                            is used need to specify the algo for sparse
%                            coding(1-@dp, 2-@dl, 3-@dch, 4-@dchperceptron)
%
%   OUTPUT:
%   out - structure containing various output variables
%       U        - the learned dictionary
%       dict_idx - the indices of points from the data matrix that are 
%                  chosen to be dictionary atoms (NA for KSVD and randomized)
%

timerVal = tic;
%rng(0);
rng('default');

P = csvread(data_path);
%%% Only for MNIST data which has been Fourier Transformed
%%%P = abs(P);
%%%
[d, n] = size(P);
[~, file_name, ~] = fileparts(data_path);
fprintf('\n');
fprintf('-------------------------------------------------------------\n');
display(['Dataset name:' file_name]);
display(['No. of points(n):' num2str(n)]);
display(['No. of dimensions(d):' num2str(d)]);

if (~isfield(param, 'epsilon'))
    closest = pdist2(P', P', 'euclidean', 'Smallest', 2);
    %epsilon = min(closest(2, :)); % Dist between closest two points
    epsilon = mean(closest(2, :)); % Avg distance between pairs of closest points
else
    epsilon = param.epsilon;
end
display(['Error tolerance (epsilon)=' num2str(epsilon)]);

if (~isfield(param, 'max_sparsity'))
    if algorithm_id > 2
        max_sparsity = d;
    else
        max_sparsity = algorithm_id;
    end
else
    max_sparsity = param.max_sparsity;
end
display(['Max sparsity (m)=' num2str(max_sparsity)]);

if (~isfield(param, 'stopping_criterion'))
    stopping_criterion = 1;
else
    stopping_criterion = param.stopping_criterion;
end

stopping_func_list = {@max, @mean};
stopping_func = stopping_func_list{stopping_criterion};
display(['Stopping criterion=' func2str(stopping_func)]);

if algorithm_id == 5 || algorithm_id == 6
    if (~isfield(param, 'dictionary_size'))
        k = 50;
    else
        k = param.dictionary_size;
    end
end

if algorithm_id == 6
    assert(isfield(param, 'sparse_algo'), 'Algo for sparse coding unspecified!');
end

algorithm_list = {@dp, @dl, @dch, @dch, @KSVD, 'randomized'};
algorithm = algorithm_list{algorithm_id};
switch algorithm_id
    case 1
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        [selected, ~, dist_array, avg_dist_array, count_inactive] = algorithm(P, epsilon, stopping_func);        
        U = P(:, selected);
        out.U = U;
        out.dict_idx = selected;
    case 2
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, epsilon, stopping_func);        
        [~, dict_idx] = ismember(U', P', 'rows');
        out.U = U;
        out.dict_idx = dict_idx;
    case 3
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, 1, epsilon, max_sparsity, stopping_func);        
        [~, dict_idx] = ismember(U', P', 'rows');
        out.U = U;
        out.dict_idx = dict_idx;
    case 4
        algorithm_name = [func2str(algorithm) 'perceptron'];
        display(['Algorithm chosen:' algorithm_name]);
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, 2, epsilon, max_sparsity, stopping_func);
        [~, dict_idx] = ismember(U', P', 'rows');
        out.U = U;
        out.dict_idx = dict_idx;
    case 5
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        param1.K = k;
        param1.numIteration = 50;
        param1.errorFlag = 1;
        param1.preserveDCAtom = 0;
        param1.errorGoal = epsilon;
        param1.InitializationMethod =  'DataElements';
        param1.displayProgress = 1;
        [U, output] = algorithm(P, param1);
        out.U = U;
    case 6
        display(['Algorithm chosen:' algorithm]);
        U = datasample(P, k, 2, 'Replace', false);
        [~, dict_idx] = ismember(U', P', 'rows');
        out.U = U;
        out.dict_idx = dict_idx;
    otherwise
        disp('Incorrect input');
        return
end

k = size(U, 2);
disp(['Size of dictionary learned:' num2str(k)]);

% Plot the distance as a function of dictionary size (only for dp, dl and
% dch)
% if algorithm_id < 5
%     % Plot and save distance (max and average) as a function of the size of 
%     % the dictionary size
%     distance_figure = figure('visible', 'off');
% 
%     subplot(2, 1, 1);
%     plot(dist_array);
%     title([file_name '\_' algorithm_name '\_distance\_epsilon=' num2str(epsilon) '\_maxsparsity=' num2str(max_sparsity) '\_stoppingcriterion=' func2str(stopping_func)]);
%     xlabel('No. of points in U');
%     hold on;
%     plot(avg_dist_array, 'green');
%     legend('Max dist', 'Avg dist');
%     hold off;
% 
%     subplot(2, 1, 2);
%     plot(count_inactive);
%     title(['Count of total inactive points in (P + U) (<= epsilon), epsilon = ' num2str(epsilon)]);
%     xlabel('No. of points in U');
%     legend('Inactive points', 'Location', 'NorthWest');
% 
%     saveas(distance_figure, ['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_distance_epsilon=' num2str(epsilon) '_maxsparsity=' num2str(max_sparsity) '_stoppingcriterion=' func2str(stopping_func) '.jpg']);
% end

% Get the sparsity level when @dch is used and plot the average distance vs
% sparsity level. Note that for @dp and @dl, sparsity level is one and two 
% respectively.
if algorithm_id == 6
    algorithm_id = param.sparse_algo;
    algorithm_name = [algorithm '-' func2str(algorithm_list{algorithm_id})];
    if algorithm_id == 4
        algorithm_name = [algorithm_name 'perceptron'];
    end
end

if algorithm_id < 5
    method_list = {@pdist2, @compute_dist_closest_line, @compute_dist_chull, @compute_dist_chull_perceptron};
    method = method_list{algorithm_id};
end
switch algorithm_id
    case 1
        sparsity_level = 1;
        memory_final_1 = k*(d-1) + n;
        D = method(U', P', 'Euclidean', 'Smallest', 1);
        final_cost = stopping_func(D);
        sparsity_coeff = (k + (n-k)*sparsity_level)/n;
        memory_final_2 = d*k + n*sparsity_coeff;
    case 2
        sparsity_level = 2;
        memory_final_1 = k*(d-3) + 3*n;
        D = method(U, P);
        final_cost = stopping_func(D);
        sparsity_coeff = (k + (n-k)*sparsity_level)/n;
        memory_final_2 = d*k + n*sparsity_coeff;
    case 5
        X = output.CoefMatrix;
        sparsity_level = k;
        sparsity_coeff = output.numCoef(end);
        max_sparsity = k;
        memory_final_1 = d*k + 2*nnz(X);
        memory_final_2 = d*k + nnz(X);
        error_matrix = P - U*X;
        final_cost = mean(sqrt(sum(error_matrix.^2)));
    otherwise
        %avg_dist_with_sparsity = zeros(1, max_sparsity);
        sparsity_level = max_sparsity;
        flag = 0;
        sparsity_count = zeros(1, max_sparsity);
        previous_dist = [];
        for j = 1:max_sparsity            
            disp(['Sparsity iteration:' num2str(j)]);
            
            D = method(U, P, j);
            %avg_dist_with_sparsity(j) = mean(D);
            if stopping_criterion == 1
                threshold = epsilon;
            else
                threshold = 2*epsilon - max(D);
            end
            threshold = max(threshold, 0);
            
            previous_dist = [previous_dist D(D <= threshold)];
            sparsity_count(j) = sum(D <= threshold);
            P = P(:, D > threshold);
            D = D(D > threshold);            
            all_dist = [previous_dist D];            
            final_cost = stopping_func(all_dist);
            if final_cost <= epsilon && flag == 0
                sparsity_level = j;
                flag = 1;
                sparsity_count(j) = sparsity_count(j) + size(D, 2);
                break;
            end
        end
        sparsity_count = sparsity_count(1:sparsity_level);
        if flag == 0
            sparsity_count(sparsity_level) = sparsity_count(sparsity_level) + size(D, 2);
        end
        disp(sparsity_count);
        disp(size(D, 2));
        disp(sum(sparsity_count));
        assert(sum(sparsity_count) == n, 'Sparse approximation not generated for all points!');
        sparsity_coeff = (sum((1:sparsity_level).*sparsity_count))/n;
        
%         avg_dist_vs_sparsity_figure = figure('visible', 'off');
%         plot(avg_dist_with_sparsity);
%         title([file_name '\_' algorithm_name '\_avgdistance\_vs\_sparsity\_epsilon=' num2str(epsilon) '\_maxsparsity=' num2str(max_sparsity) '\_stoppingcriterion=' func2str(stopping_func)]);
%         xlabel('Sparsity level');
%         ylabel('Avg. distance (error)');
%         refline(0, epsilon);
%         saveas(avg_dist_vs_sparsity_figure, ['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_avgdistance_vs_sparsity_epsilon=' num2str(epsilon) '_maxsparsity=' num2str(max_sparsity) '_stoppingcriterion=' func2str(stopping_func) '.jpg']);
%         
        if algorithm_id == 3
            %memory_final = k*d + (2*sparsity_level-1)*(n-k);
            memory_final_1 = k*d + sparsity_count(1)-k + sum((2*(2:sparsity_level)-1).*sparsity_count(2:end));
            memory_final_2 = d*k + n*sparsity_coeff;
        end
        if algorithm_id == 4
            %memory_final = k*d + sparsity_level*(n-k);
            memory_final_1 = k*d + sparsity_count(1)-k + sum((2*(2:sparsity_level)-1).*sparsity_count(2:end));
            memory_final_2 = d*k + n*sparsity_coeff;
        end
end        

% Compute various performance metrics, print them, and store the results
memory_initial = n*d;
compression_ratio_1 = memory_final_1/memory_initial;
compression_ratio_2 = memory_final_2/memory_initial;
results = [n, d, epsilon, k, memory_initial, memory_final_1, memory_final_2, compression_ratio_1, compression_ratio_2, sparsity_level, sparsity_coeff, final_cost, max_sparsity];
string = sprintf('%0.5f,', results);
string = [string func2str(stopping_func)];
% if exist(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_performance_metrics_epsilon=' num2str(epsilon) '_maxsparsity=' num2str(max_sparsity) '_stoppingcriterion=' func2str(stopping_func) '.csv'], 'file') == 0
%     header = 'No. of points(n),No. of dimensions(d),Error tolerance(epsilon),Dictionary size(k),Initial memory,Final memory,Compression ratio,Actual max sparsity,Sparsity coeff,Final cost,Max allowed sparsity(m),Stopping criterion\n';
%     fid = fopen(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_performance_metrics_epsilon=' num2str(epsilon) '_maxsparsity=' num2str(max_sparsity) '_stoppingcriterion=' func2str(stopping_func) '.csv'], 'w');
%     fprintf(fid, header);
%     fclose(fid);
% end
% fid = fopen(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_performance_metrics_epsilon=' num2str(epsilon) '_maxsparsity=' num2str(max_sparsity) '_stoppingcriterion=' func2str(stopping_func) '.csv'], 'a');
% fprintf(fid, string);
% fclose(fid);

string = [file_name ',' algorithm_name ',' string '\n'];
fid = fopen('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\results14.csv', 'a');
fprintf(fid, string);
fclose(fid);    

disp(['Sparsity level:' num2str(sparsity_level)]);
disp(['Sparsity coefficient:' num2str(sparsity_coeff)]);
disp(['Final cost:' num2str(final_cost)]);
disp(['Compression ratio1:' num2str(compression_ratio_1)]);
disp(['Compression ratio2:' num2str(compression_ratio_2)]);



% % ------- Plotting the average error as a function of memory used -----
% k = size(U, 2);
% memory = (1:k)*(d-2*sparsity_level + 1) + (2*sparsity_level - 1)*n;
% Memory_plot = figure('visible', 'off');
% plot(memory, dist_array);
% title(['Error vs Memory Chull-wise Perceptron, epsilon = ' num2str(epsilon)]);
% xlabel('Memory');
% hold on;
% plot(memory, avg_dist_array, 'green');
% legend('Max dist', 'Avg dist');
% hold off;
% saveas(Memory_plot, ['Output-' file_name '\Perceptron\chull-wise-memory' num2str(iter) '_epsilon=' num2str(epsilon) '.jpg']);
%
%     % Evaluate X using matrix multiplication
%     P = Q; % Since P was changed in computing U, we need to reaasign
%     c = 0.1;
%     U = U + c*eye(size(U));
%     X = U\P;
%
%     % Cost (Frobenius norm) ||P - U*X||^2
%     C = norm((P - U*X), 'fro');
%
%     % Sparsity
%     nonzero_count = sum(sum(abs(X) > 1e-5));
%     sparsity_coeff = nonzero_count/numel(X);
%
%     % % Output
%     % fprintf('No. of points(n): %d\t', n);
%     % fprintf('No. of dimensions(d): %d\n', d);
%     % fprintf('No. of clusters(k): %d\n', size(U, 2));
%     % fprintf('Cost ||P - UX||: %3.3f\n', C);
%     % fprintf('Sparsity: %3.3f\n', sparsity_coeff);
%
%     % Output - saved in csv file
%     output = [n, d, size(U, 2), C, sparsity_coeff, epsilon, iterations];
%     % ------ CHANGE HERE(specify directory)-----------
%     dlmwrite('Output-wdbc\Max_avg_inactive\Chull-wise\results_chull.csv', output, '-append');
%end
toc(timerVal);
end

