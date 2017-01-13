function [] = normalize_dataset(data_path)
% NORMALIZE_DATASET Preprocess the dataset to have all feature values 
% between 0 and 1.
%
%   INPUT:
%   data_path - full path to the chosen dataset
%
%   OUTPUT:
%   No variables returned. Saves the modified dataset at the default
%   location.
%
P = csvread(data_path);
[~, file_name, ~] = fileparts(data_path);

[~, d] = size(P);
for i = 1:d
    if max(P(:, i)) == min(P(:, i))
        P(:, i) = 1;
    else
        P(:, i) = (P(:, i) - min(P(:, i)))/(max(P(:, i)) - min(P(:, i)));
    end
end
P = P';

dlmwrite(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\' file_name '-norm.csv'], P);
end

