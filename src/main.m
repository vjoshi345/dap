function [] = main(data_path, algorithm_id, epsilon)
% MAIN Function to perform dictionary learning on a dataset with a user
% specified algorithm.
%
%   INPUT:
%   data_path    - full path to the dataset stored as a .csv or .mat file
%   algorithm_id - algorithm to perform dictionary learning (specified as a
%                  number)
%                  1 - distance from points
%                  2 - distance from line segments
%                  3 - distance from convex hull
%                  4 - distance from convex hull (perceptron method)
%   epsilon      - optional argument specifying the error tolerance
%
%   OUTPUT:
%   No variables returned. Summary statistics are plotted and saved.
%
rng(0);

P = csvread(data_path);
P = P';
[~, file_name, ~] = fileparts(data_path);
display(['Dataset name:' file_name]);
display(['No. of points(n):' num2str(size(P, 2))]);
display(['No. of dimensions(d):' num2str(size(P, 1))]);

if nargin < 3
    % Choosing epsilon
    closest = pdist2(P', P', 'euclidean', 'Smallest', 2);
    %epsilon = min(closest(2, :)); % Dist between closest two points
    epsilon = mean(closest(2, :)); % Avg distance between pairs of closest points
end
display(['Error tolerance (epsilon)=' num2str(epsilon)]);

algorithm_list = {@dp, @dl, @dch, @dch};
algorithm = algorithm_list{algorithm_id};

switch algorithm_id
    case 1 | 2
        [selected, sparse_code, dist_array, avg_dist_array, count_inactive] = algorithm(P, epsilon);
        algorithm_name = func2str(algorithm);
    case 3
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, 1, epsilon);
        algorithm_name = func2str(algorithm);
    case 4
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, 2, epsilon);
        algorithm_name = [func2str(algorithm) 'perceptron'];
    otherwise
        disp('Incorrect input');
        return        
end
display(['Algorithm chosen:' algorithm_name]);

% Plot and save distance as a function of the size of subset chosen
distance_figure = figure('visible', 'off');
%distance_figure = figure();

subplot(2, 1, 1);
plot(dist_array);
title([file_name '\_' algorithm_name '\_distance\_epsilon=' num2str(epsilon)]);
xlabel('No. of points in U');
hold on;
plot(avg_dist_array, 'green');
legend('Max dist', 'Avg dist');
hold off;

subplot(2, 1, 2);
plot(count_inactive);
title(['Count of total inactive points in (P + U) (<= epsilon), epsilon = ' num2str(epsilon)]);
xlabel('No. of points in U');
legend('Inactive points', 'Location', 'NorthWest');

saveas(distance_figure, ['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_distance_epsilon=' num2str(epsilon) '.jpg']);


% % -------- Approximating points in P with various sparsity levels -----
% dist_with_sparsity = zeros(1, d);
% sparsity_level = Inf;
% for j = 1:d
%     %D = distance_chull(U, P, j);
%     D = distance_chull_perceptron(U, P, j);
%     dist_with_sparsity(j) = mean(D);
%     if mean(D) <= epsilon && sparsity_level == Inf
%         sparsity_level = j;
%     end
% end
% 
% Avg_cost = figure('visible', 'off');
% 
% plot(dist_with_sparsity);
% title(['Avg. distance of points in P at different sparsity levels, epsilon = ' num2str(epsilon)]);
% xlabel('Sparsity level');
% ylabel('Avg. distance');
% refline(0, epsilon);
% 
% % ------ CHANGE HERE(specify directory)-----------
% saveas(Avg_cost, ['Output-' file_name '\Perceptron\Cost_with_sparsity' num2str(iter) '_epsilon=' num2str(epsilon) '.jpg']);
% 
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
% % ------ Computing sparsity coeff and cost and storing results --------
% C = dist_with_sparsity(sparsity_level);
% k = size(U, 2);
% sparsity_coeff = (k + (n-k)*sparsity_level)/(n*k);
% output = [n, d, epsilon, k, sparsity_level, sparsity_coeff, C];
% 
% memory_initial = n*d;
% memory_final = memory(end);
% compression_ratio = memory_final/memory_initial;
% memory_output = [n, d, epsilon, k, memory_initial, memory_final, compression_ratio, sparsity_level, sparsity_coeff, C];
% % ------ CHANGE HERE(specify directory)-----------
% dlmwrite(['Output-' file_name '\Perceptron\savings.csv'], output, '-append');
% dlmwrite(['Output-' file_name '\Perceptron\memory.csv'], memory_output, '-append');
% 
% display(iter);

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




end

