function [U, dist_array, avg_dist_array, count_inactive] = dp(P, epsilon)
% DP Learns a dictionary for an input dataset using the distance from 
% points algorithm.
%
%   INPUT:
%   P       - the input dataset as a matrix with columns as datapoints 
%             and rows as dimensions
%   epsilon - error tolerance for each datapoint
%   
%   OUTPUT:
%   U              - the dictionary learned by the algorithm as a matrix
%   dist_array     - distance of farthest point from dictionary atoms at 
%                    each iteration of the algorithm
%   avg_dist_array - average distance of points from dictionary atoms at 
%                    each iteration of the algorithm
%   count_inactive - number of points within distance epsilon from 
%                    dictionary atoms at each iteration of the algorithm
%
% P = U*X
% P = d*n matrix, U = d*k matrix, X = k*n matrix
% 

%tic;

% for iter = 1:1
%     clearvars -except iter;
%     % Setting Constants
%     %epsilon = iter;
% 
%     % Importing data and converting to the matrix form
%     file_name = '10line-data'; % ----- CHANGE HERE (specify file)-------------
%     P = csvread([file_name '-mod.csv']); 
%     P = P'; 
%     %P = P(:, 1:50);
    [d, n] = size(P);
    
    if nargin < 2
        % Choosing epsilon
        closest = pdist2(P', P', 'euclidean', 'Smallest', 2);
        %epsilon = min(closest(2, :)); % Dist between closest two points
        epsilon = mean(closest(2, :)); % Avg dist between closest pair of points
    end
    
    
    %Q = P; % We will manipuate P and keep a copy of it in Q for later use

    %%%%%%%%%%%%%%%%%% ALGORITHM-1 Point-wise max distance %%%%%%%%%%%%%%%%%%%%
    % Computing columns of U 
    r = randi([1, n], 1, 1);
    U = zeros(d, n);
    U(:, 1) = P(:, r);
    %P = P(:, [(1:r-1), (r+1:end)]);
    P(:, r) = [];

    dist_array = zeros(1, n);
    [D, ~] = pdist2((U(:, 1))', P', 'euclidean', 'Smallest', 1);
    [max_dist, max_index] = max(D);
    dist_array(1) = max_dist;

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
                [D, ~] = pdist2((U(:, 1:i))', P', 'euclidean', 'Smallest', 1);
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
%     Point_wise = figure('visible', 'off');
% 
%     subplot(2, 1, 1);
%     plot(dist_array);
%     title(['Point-wise distance, epsilon = ' num2str(epsilon)]);
%     xlabel('No. of points in U');
%     hold on;
%     plot(avg_dist_array, 'green');
%     legend('Max dist', 'Avg dist');
%     hold off;
% 
%     subplot(2, 1, 2);
%     plot(count_inactive);
%     title(['Count of inactive points in P (<= epsilon), epsilon = ' num2str(epsilon)]);
%     xlabel('No. of points in U');
%     legend('Inactive points');
%     
%     % ------- Plotting the average error as a function of memory used -----
%     k = size(U, 2);
%     memory = (1:k)*(d-1) + n;
%     Memory_plot = figure('visible', 'off');
%     plot(memory, dist_array);
%     title(['Error vs Memory Point-wise, epsilon = ' num2str(epsilon)]);
%     xlabel('Memory');
%     hold on;
%     plot(memory, avg_dist_array, 'green');
%     legend('Max dist', 'Avg dist');
%     hold off;
%     
%     % ------ CHANGE HERE(specify directory)-----------
%     saveas(Point_wise, ['Output-' file_name '\Point-wise\point_wise' num2str(iter) '_epsilon=' num2str(epsilon) '.jpg']); 
%     saveas(Memory_plot, ['Output-' file_name '\Point-wise\memory' num2str(iter) '_epsilon=' num2str(epsilon) '.jpg']); 
%     
%     % ------ Computing sparsity coeff and cost and storing results --------
%     C = avg_dist_array(end);
%     k = size(U, 2);
%     sparsity_level = 1; % Since this is the point-wise algorithm
%     sparsity_coeff = (k + (n-k)*sparsity_level)/(n*k);
%     output = [n, d, epsilon, k, sparsity_level, sparsity_coeff, C];
%     
%     memory_initial = n*d;
%     memory_final = memory(end);
%     compression_ratio = memory_final/memory_initial;
%     memory_output = [n, d, epsilon, k, memory_initial, memory_final, compression_ratio, sparsity_level, sparsity_coeff, C];
%     % ------ CHANGE HERE(specify directory)-----------
%     dlmwrite(['Output-' file_name '\Point-wise\savings.csv'], output, '-append');
%     
%     dlmwrite(['Output-' file_name '\Point-wise\memory.csv'], memory_output, '-append');
%     
%     display(iter);
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
%     output = [n, d, r, size(U, 2), C, sparsity_coeff, epsilon];
%     % ------ CHANGE HERE(specify directory)-----------
%     dlmwrite('Output-iono\Max_avg_inactive\Point-wise\results_point.csv', output, '-append');
end


%toc;


