function [U, dist_array, avg_dist_array, count_inactive] = dch(P, epsilon, iterations)
% DCH Learns a dictionary for an input dataset using the distance from 
% convex-hull algorithm.
%
%   INPUT:
%   P          - the input dataset as a matrix with columns as datapoints 
%                and rows as dimensions.
%   epsilon    - error tolerance for each datapoint.
%   iterations - number of iterations to compute distance to convex hull.
%   
%   OUTPUT:
%   U              - the dictionary learned by the algorithm as a matrix.
%   dist_array     - distance of farthest point from the convex-hull of the
%                    dictionary at each iteration of the algorithm.
%   avg_dist_array - average distance of points from the convex-hull of the
%                    dictionary at each iteration of the algorithm.
%   count_inactive - number of points within distance epsilon from the
%                    convex-hull of the dictionary at each iteration.
%   

tic;
%for iter = 1:1
    %clearvars -except iter;
    % Set seed (do this if you want to choose the same initial point) -----
    %rng(0);
    
    % Setting Constants
    %epsilon = 1;
    %iterations = 10;

    % Importing data and converting to the matrix form
    % ----- CHANGE HERE (specify file)---------------
%     file_name = 'CTG';
%     P = csvread([file_name '-mod.csv']);
%     P = P'; 
%     %P = P(:, 1:150);
    [d, n] = size(P);
    
    if nargin < 3
        iterations = d;
    end
    if nargin < 2
        % Choosing epsilon
        closest = pdist2(P', P', 'euclidean', 'Smallest', 2);
        %epsilon = min(closest(2, :)); % Dist between closest two points
        epsilon = mean(closest(2, :)); % Avg distance between pairs of closest points
    
        %epsilon = (iter+1)*epsilon/2; % Multiplying epsilon to get results for different values
    end
    
    %Q = P; % We will manipuate P and keep a copy of it in Q for later use


    %--- DICTIONARY LEARNING - ALGORITHM-1 Distance from convex hull ------
    % Computing columns of U (greedy algorithm)
    r = randi([1, n], 1, 1);
    U = zeros(d, n);
    U(:, 1) = P(:, r);
    P(:, r) = [];

    D = distance_chull(U(:, 1), P, iterations);
    %D = distance_chull_perceptron(U(:, 1), P, iterations);
    [max_dist, max_index] = max(D);
    dist_array = zeros(1, n);
    dist_array(1) = max_dist;
    
    %display(D);
    %display(max_dist);
    
    avg_dist_array = zeros(1, n);
    avg_dist_array(1) = mean(D);

    count_inactive = zeros(1, n);
    count_inactive(1) = sum(D <= epsilon) + 1;

    flag = 0;
    for i = 2:n
        if max_dist <= epsilon
            flag = 1;
            break
        else
            U(:, i) = P(:, max_index);
            P(:, max_index) = [];
            if i == n
                dist_array(i) = 0;
                avg_dist_array(i) = 0;
                count_inactive(i) = n;
            else
                D = distance_chull(U(:, 1:i), P, iterations);
                %D = distance_chull_perceptron(U(:, 1:i), P, iterations);
                [max_dist, max_index] = max(D);
                dist_array(i) = max_dist;
                avg_dist_array(i) = mean(D);
                count_inactive(i) = sum(D <= epsilon) + i;
                
            end
        end
        fprintf('End of iteration:%d\n', i);
    end

    if flag == 1
        U = U(:, 1:(i-1));
        dist_array  = dist_array(1:(i-1));
        avg_dist_array  = avg_dist_array(1:(i-1));
        count_inactive = count_inactive(1:(i-1));
    else
        U = U(:, 1:i);
        dist_array = dist_array(1:i);
        avg_dist_array  = avg_dist_array(1:i);
        count_inactive = count_inactive(1:i);
    end

%     %%% Plotting the distance as a function of the size of subset chosen
%     Chull_wise = figure('visible', 'off');
% 
%     subplot(2, 1, 1);
%     plot(dist_array);
%     title(['Chull-wise distance, epsilon = ' num2str(epsilon)]);
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
%     
%         
%     % ------ CHANGE HERE(specify directory)-----------
%     saveas(Chull_wise, ['Output-' file_name '\Chull-ista\dict_learn' num2str(iter) '_epsilon=' num2str(epsilon) '.jpg']);
%     
%     % -------------------------- SPARSE CODING ----------------------------
%     P = Q; % reassigning Q to P to get back the original matrix
%     k = size(U, 2);
%  
%     %------------- ISTA/LASSO algorithm for sparse coding -----------------------
%     error = zeros(1, 5);
%     sparsity = zeros(1, 5);
%     cost_func = zeros(1, 5);
%     sparsity_coeff = zeros(1, 5);
%     tau = zeros(1, 5);
%     
%     for j = 1:5
%         tau(j) = 0.5*j;
%         X = zeros(k, n);
%     
% %         % ISTA for sparse coding
% %         for i = 1:n
% %             [x, ~] = TwIST(P(:, i), U, tau(j), 'Phi', @(n) sum(abs(n)));
% %             X(:, i) = x;
% %         end
%         
%         % LASSO for sparse coding
%         for i = 1:n
%             [x, fit_info] = lasso(U, P(:, i), 'Lambda', tau(j));    
%             X(:, i) = x';
%         end
%                
%         error(j) = norm(P - U*X, 'fro');
%         sparsity(j) = tau(j)*sum(sum(abs(X), 1)); % Using the L1-norm as a sparsity measure
%         cost_func(j) = error(j) + sparsity(j);
%         
%         sparsity_coeff(j) = nnz(X)/(k*n);
%         
%         memory_initial = d*n;
%         memory_final = d*k + nnz(X);
%         compression_ratio = memory_final/memory_initial;
%         
%         results = [n, d, epsilon, tau(j), k, error(j), sparsity_coeff(j), sparsity(j), cost_func(j), memory_initial, memory_final, compression_ratio];
%         dlmwrite(['Output-' file_name '\Chull-lasso\results.csv'], results, '-append');
%         
%         display(['iter:' j]);
%     end
%     
%     Chull_wise = figure('visible', 'off');
%     
%     subplot(2, 1, 1);
%     plot(tau, error);
%     title(['Reconstruction Error Chull LASSO, epsilon = ' num2str(epsilon)]);
%     xlabel('Regularization parameter (tau)');
%     
%     subplot(2, 1, 2);
%     plot(tau, sparsity_coeff, 'green');
%     title(['Sparsity Chull LASSO, epsilon = ' num2str(epsilon)]);
%     xlabel('Regularization parameter (tau)');
%     
%     saveas(Chull_wise, ['Output-' file_name '\Chull-lasso\sparse_coding' num2str(iter) '_epsilon=' num2str(epsilon) '.jpg']);
    
%     % -------- Approximating points in P with various sparsity levels -----
%     dist_with_sparsity = zeros(1, d);
%     sparsity_level = Inf;
%     for j = 1:d
%         %D = distance_chull(U, P, j);
%         D = distance_chull_perceptron(U, P, j);
%         dist_with_sparsity(j) = mean(D);
%         if mean(D) <= epsilon && sparsity_level == Inf
%             sparsity_level = j;
%         end
%     end
%     
%     Avg_cost = figure('visible', 'off');
% 
%     plot(dist_with_sparsity);
%     title(['Avg. distance of points in P at different sparsity levels, epsilon = ' num2str(epsilon)]);
%     xlabel('Sparsity level');
%     ylabel('Avg. distance');
%     refline(0, epsilon);
%     
%     % ------ CHANGE HERE(specify directory)-----------
%     saveas(Avg_cost, ['Output-' file_name '\Perceptron\Cost_with_sparsity' num2str(iter) '_epsilon=' num2str(epsilon) '.jpg']);
%     
%     % ------- Plotting the average error as a function of memory used -----
%     k = size(U, 2);
%     memory = (1:k)*(d-2*sparsity_level + 1) + (2*sparsity_level - 1)*n;
%     Memory_plot = figure('visible', 'off');
%     plot(memory, dist_array);
%     title(['Error vs Memory Chull-wise Perceptron, epsilon = ' num2str(epsilon)]);
%     xlabel('Memory');
%     hold on;
%     plot(memory, avg_dist_array, 'green');
%     legend('Max dist', 'Avg dist');
%     hold off;
%     saveas(Memory_plot, ['Output-' file_name '\Perceptron\chull-wise-memory' num2str(iter) '_epsilon=' num2str(epsilon) '.jpg']); 
%     
%     % ------ Computing sparsity coeff and cost and storing results --------
%     C = dist_with_sparsity(sparsity_level);
%     k = size(U, 2);
%     sparsity_coeff = (k + (n-k)*sparsity_level)/(n*k);
%     output = [n, d, epsilon, k, sparsity_level, sparsity_coeff, C];
%     
%     memory_initial = n*d;
%     memory_final = memory(end);
%     compression_ratio = memory_final/memory_initial;
%     memory_output = [n, d, epsilon, k, memory_initial, memory_final, compression_ratio, sparsity_level, sparsity_coeff, C];
%     % ------ CHANGE HERE(specify directory)-----------
%     dlmwrite(['Output-' file_name '\Perceptron\savings.csv'], output, '-append');
%     dlmwrite(['Output-' file_name '\Perceptron\memory.csv'], memory_output, '-append');
    
    %display(iter);
toc;    
end

%toc;

%end

