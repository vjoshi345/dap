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
[d, n] = size(P);
[~, file_name, ~] = fileparts(data_path);
display(['Dataset name:' file_name]);
display(['No. of points(n):' num2str(n)]);
display(['No. of dimensions(d):' num2str(d)]);

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
    case 1
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        [selected, ~, dist_array, avg_dist_array, count_inactive] = algorithm(P, epsilon);        
        U = P(:, selected);
    case 2
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, epsilon);        
    case 3
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, 1, epsilon);        
    case 4
        algorithm_name = [func2str(algorithm) 'perceptron'];
        display(['Algorithm chosen:' algorithm_name]);
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, 2, epsilon);        
    otherwise
        disp('Incorrect input');
        return
end

k = size(U, 2);
disp(['Size of dictionary leanred:' num2str(k)]);

% Plot and save distance (max and average) as a function of the size of 
% the dictionary size
distance_figure = figure('visible', 'off');

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

% Get the sparsity level when @dch is used and plot the average distance vs
% sparsity level. Note that for @dp and @dl, sparsity level is one and two 
% respectively.
method_list = {[], [], @compute_dist_chull, @compute_dist_chull_perceptron};
method = method_list{algorithm_id};
switch algorithm_id
    case 1
        sparsity_level = 1;
        memory_final = k*(d-1) + n;
    case 2
        sparsity_level = 2;
        memory_final = k*(d-3) + 3*n;
    otherwise
        avg_dist_with_sparsity = zeros(1, d);
        sparsity_level = Inf;
        for j = 1:d
            disp(['Sparsity iteration:' num2str(j)]);
            D = method(U, P, j);
            avg_dist_with_sparsity(j) = mean(D);
            if mean(D) <= epsilon && sparsity_level == Inf
                sparsity_level = j;
            end
        end
        avg_dist_vs_sparsity_figure = figure('visible', 'off');
        plot(avg_dist_with_sparsity);
        title([file_name '\_' algorithm_name '\_avgdistance\_vs\_sparsity\_epsilon=' num2str(epsilon)]);
        xlabel('Sparsity level');
        ylabel('Avg. distance (error)');
        refline(0, epsilon);
        saveas(avg_dist_vs_sparsity_figure, ['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_avgdistance_vs_sparsity_epsilon=' num2str(epsilon) '.jpg']);
        if algorithm_id == 3
            memory_final = k*d + (2*sparsity_level-1)*(n-k);
        end
        if algorithm_id == 4
            memory_final = k*d + sparsity_level*(n-k);
        end
end        

% Compute various performance metrics, print them, and store the results
final_cost = avg_dist_with_sparsity(sparsity_level);
sparsity_coeff = (k + (n-k)*sparsity_level)/(n*k);
memory_initial = n*d;
compression_ratio = memory_final/memory_initial;
results = [n, d, epsilon, k, memory_initial, memory_final, compression_ratio, sparsity_level, sparsity_coeff, final_cost];
if exist(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_performance_metrics_epsilon=' num2str(epsilon) '.csv'], 'file') == 0
    header = 'No. of points(n),No. of dimensions(d),Error tolerance(epsilon),Dictionary size(k),Initial memory,Final memory,Compression ratio,Sparsity level,Sparsity coeff,Final cost\n';
    fid = fopen(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_performance_metrics_epsilon=' num2str(epsilon) '.csv'], 'w');
    fprintf(fid, header);
    fclose(fid);
end
dlmwrite(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_performance_metrics_epsilon=' num2str(epsilon) '.csv'], results, '-append');

disp(['Sparsity level:' num2str(sparsity_level)]);
disp(['Sparsity coefficient:' num2str(sparsity_coeff)]);
disp(['Final cost:' num2str(final_cost)]);
disp(['Compression ratio:' num2str(compression_ratio)]);



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

end

