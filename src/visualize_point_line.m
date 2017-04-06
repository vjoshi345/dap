function [] = visualize_point_line(data_path, stopping_func)
rng('default');
Y = csvread(data_path);

closest = pdist2(Y', Y', 'euclidean', 'Smallest', 2);
epsilon = mean(closest(2, :)); % Avg distance between pairs of closest points

[~, sparse_code, ~, ~, ~] = dp(Y, epsilon, stopping_func);        
uv = unique(sparse_code);
counts = zeros(length(uv), 1);
for i = 1:length(uv)
    counts(i) = sum(sparse_code == uv(i));
end

[sorted_counts, index] = sort(counts, 1, 'descend');
% atom_idx = uv(index(1:10));
% atom_idx = cellstr(num2str(atom_idx(:)));
% c = categorical(atom_idx);
bar(sorted_counts)

% sparse_code = categorical(sparse_code);
% summary(sparse_code)
% sort(sparse_code)

end