function [] = main(data_path, algorithm_id, stopping_criterion, max_sparsity, epsilon)
% MAIN Function to perform dictionary learning on a dataset with a user
% specified algorithm.
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
%   max_sparsity       - optional argument specifying the no. of iterations
%                        to compute the distance to convex hull for the dch 
%                        algorithms (default = #dimensions of input data)
%   stopping_criterion - optional argument specifying either the max(1) or
%                        mean(2) as the stopping criteria (default = max)
%   epsilon            - optional argument specifying the error tolerance
%
%   OUTPUT:
%   No variables returned. Summary statistics are plotted and saved.
%
timerVal = tic;
rng(0);

P = csvread(data_path);
[d, n] = size(P);
[~, file_name, ~] = fileparts(data_path);
display(['Dataset name:' file_name]);
display(['No. of points(n):' num2str(n)]);
display(['No. of dimensions(d):' num2str(d)]);

if nargin < 5
    % Choosing epsilon
    closest = pdist2(P', P', 'euclidean', 'Smallest', 2);
    %epsilon = min(closest(2, :)); % Dist between closest two points
    epsilon = mean(closest(2, :)); % Avg distance between pairs of closest points
end
display(['Error tolerance (epsilon)=' num2str(epsilon)]);

if nargin < 4 
    if algorithm_id > 2
        max_sparsity = d;    
    else
        max_sparsity = algorithm_id;
    end
end
display(['Max sparsity (m)=' num2str(max_sparsity)]);

if nargin < 3
    stopping_criterion = 1;
end
stopping_func_list = {@max, @mean};
stopping_func = stopping_func_list{stopping_criterion};
display(['Stopping criterion=' func2str(stopping_func)]);

algorithm_list = {@dp, @dl, @dch, @dch};
algorithm = algorithm_list{algorithm_id};

switch algorithm_id
    case 1
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        [selected, ~, dist_array, avg_dist_array, count_inactive] = algorithm(P, epsilon, stopping_func);        
        U = P(:, selected);
    case 2
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, epsilon, stopping_func);        
    case 3
        algorithm_name = func2str(algorithm);
        display(['Algorithm chosen:' algorithm_name]);
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, 1, epsilon, max_sparsity, stopping_func);        
    case 4
        algorithm_name = [func2str(algorithm) 'perceptron'];
        display(['Algorithm chosen:' algorithm_name]);
        [U, dist_array, avg_dist_array, count_inactive] = algorithm(P, 2, epsilon, max_sparsity, stopping_func);        
    otherwise
        disp('Incorrect input');
        return
end

k = size(U, 2);
disp(['Size of dictionary learned:' num2str(k)]);

% Plot and save distance (max and average) as a function of the size of 
% the dictionary size
distance_figure = figure('visible', 'off');

subplot(2, 1, 1);
plot(dist_array);
title([file_name '\_' algorithm_name '\_distance\_epsilon=' num2str(epsilon) '\_maxsparsity=' num2str(max_sparsity) '\_stoppingcriterion=' func2str(stopping_func)]);
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

saveas(distance_figure, ['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_distance_epsilon=' num2str(epsilon) '_maxsparsity=' num2str(max_sparsity) '_stoppingcriterion=' func2str(stopping_func) '.jpg']);

% Get the sparsity level when @dch is used and plot the average distance vs
% sparsity level. Note that for @dp and @dl, sparsity level is one and two 
% respectively.
method_list = {@pdist2, @compute_dist_closest_line, @compute_dist_chull, @compute_dist_chull_perceptron};
method = method_list{algorithm_id};
switch algorithm_id
    case 1
        sparsity_level = 1;
        memory_final = k*(d-1) + n;
        D = method(U', P', 'Euclidean', 'Smallest', 1);
        final_cost = stopping_func(D);
        sparsity_coeff = (k + (n-k)*sparsity_level)/n;
    case 2
        sparsity_level = 2;
        memory_final = k*(d-3) + 3*n;
        D = method(U, P);
        final_cost = stopping_func(D);
        sparsity_coeff = (k + (n-k)*sparsity_level)/n;
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
            memory_final = k*d + sparsity_count(1)-k + sum((2*(2:sparsity_level)-1).*sparsity_count(2:end));
        end
        if algorithm_id == 4
            %memory_final = k*d + sparsity_level*(n-k);
            memory_final = k*d + sparsity_count(1)-k + sum((2*(2:sparsity_level)-1).*sparsity_count(2:end));
        end
end        

% Compute various performance metrics, print them, and store the results
memory_initial = n*d;
compression_ratio = memory_final/memory_initial;
results = [n, d, epsilon, k, memory_initial, memory_final, compression_ratio, sparsity_level, sparsity_coeff, final_cost, max_sparsity];
string = sprintf('%0.5f,', results);
string = [string func2str(stopping_func)];
if exist(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_performance_metrics_epsilon=' num2str(epsilon) '_maxsparsity=' num2str(max_sparsity) '_stoppingcriterion=' func2str(stopping_func) '.csv'], 'file') == 0
    header = 'No. of points(n),No. of dimensions(d),Error tolerance(epsilon),Dictionary size(k),Initial memory,Final memory,Compression ratio,Sparsity level,Sparsity coeff,Final cost,Max sparsity(m),Stopping criterion\n';
    fid = fopen(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_performance_metrics_epsilon=' num2str(epsilon) '_maxsparsity=' num2str(max_sparsity) '_stoppingcriterion=' func2str(stopping_func) '.csv'], 'w');
    fprintf(fid, header);
    fclose(fid);
end
fid = fopen(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\' file_name '_' algorithm_name '_performance_metrics_epsilon=' num2str(epsilon) '_maxsparsity=' num2str(max_sparsity) '_stoppingcriterion=' func2str(stopping_func) '.csv'], 'a');
fprintf(fid, string);
fclose(fid);

string = [file_name ',' algorithm_name ',' string '\n'];
fid = fopen('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\results4.csv', 'a');
fprintf(fid, string);
fclose(fid);    

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
toc(timerVal);
end

