% P = U*X
% P = d*n matrix, U = d*k matrix, X = k*n matrix
% This program plots the max distance from the subset of points (U) chosen
% using three different algorithms as a function of the size of the subset.

tic;
for iter = 1:1
    clearvars -except iter;
    % Setting Constants
    epsilon = 1;

    % Importing data and converting to the matrix form
    P = csvread('ionosphere_mod.csv');
    P = P'; 
    [d, n] = size(P);

    Q = P; % We will manipuate P and keep a copy of it in Q for later use

    %%%%%%%%%%%%%%%%%% ALGORITHM-2 Line-wise max distance %%%%%%%%%%%%%%%%%%%%
    % Computing columns of U 
    r = randi([1, n], 1, 1);
    U = zeros(d, n);
    U(:, 1) = P(:, r);
    P(:, r) = [];

    dist_array = zeros(1, n);
    %dist_array(1) = pdist2(P', (U(:, 1))', 'euclidean', 'Largest', 1);
    %dist_array(1) = point_to_line(P, U(:, 1), U(:, 1));
    %[max_dist, max_index] = distance_line(U(:, 1), P);
    D = distance_line(U(:, 1), P); % Set of distances from each point in P to 
                                   % to the closest line-seg in U
    [max_dist, max_index] = max(D);                               
    dist_array(1) = max_dist;

    avg_dist_array = zeros(1, n);
    avg_dist_array(1) = mean(D);

    count_inactive = zeros(1, n);
    count_inactive(1) = sum(D <= epsilon);

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
                D = distance_line(U(:, 1:i), P);
                %[max_dist, max_index] = distance_line(U(:, 1:i), P);
                [max_dist, max_index] = max(D);
                dist_array(i) = max_dist;
                avg_dist_array(i) = mean(D);
                count_inactive(i) = sum(D <= epsilon);
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

    %%% Plotting the distance as a function of the size of subset chosen
    Line_wise = figure('visible', 'off');

    subplot(2, 1, 1);
    plot(dist_array);
    title(['Line-wise distance, epsilon = ' num2str(epsilon)]);
    xlabel('No. of points in U');
    hold on;
    plot(avg_dist_array, 'green');
    legend('Max dist', 'Avg dist');
    hold off;

    subplot(2, 1, 2);
    plot(count_inactive);
    title(['Count of inactive points in P (<= epsilon), epsilon = ' num2str(epsilon)]);
    xlabel('No. of points in U');
    legend('Inactive points', 'Location', 'NorthWest');
    display(iter);
    saveas(Line_wise, ['Output\Max_avg_inactive\Line-wise' num2str(iter) '_epsilon=' num2str(epsilon) '.jpg']);

    % Evaluate X using matrix multiplication
    P = Q; % Since P was changed in computing U, we need to reaasign
    c = 0.1;
    U = U + c*eye(size(U));
    X = U\P;

    % Cost (Frobenius norm) ||P - U*X||^2 
    C = norm((P - U*X), 'fro');

    % Sparsity
    nonzero_count = sum(sum(abs(X) > 1e-5));
    sparsity_coeff = nonzero_count/numel(X);

    % % Output
    % fprintf('No. of points(n): %d\t', n);
    % fprintf('No. of dimensions(d): %d\n', d);
    % fprintf('No. of clusters(k): %d\n', size(U, 2));
    % fprintf('Cost ||P - UX||: %3.3f\n', C);
    % fprintf('Sparsity: %3.3f\n', sparsity_coeff);

    % Output - saved in csv file
    output = [n, d, size(U, 2), C, sparsity_coeff, epsilon];
    dlmwrite('Output\Max_avg_inactive\results_line.csv', output, '-append');
end

toc;