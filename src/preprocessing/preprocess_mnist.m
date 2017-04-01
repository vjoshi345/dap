function [] = preprocess_mnist(dataset_size, iter)
% Script to preprocess and store the MNIST images in the data folder
% Stores a random subset of the images
rng(0);

% Read MNIST image data and get a random subset of it
oldFolder = cd('C:\CMU\CMU-Spring-2016\DAP\Dataset\MNIST database');
images = loadMNISTImages('train-images.idx3-ubyte');
labels = loadMNISTLabels('train-labels.idx1-ubyte');
[~, n] = size(images);
index = datasample(1:n, dataset_size, 'Replace', false);
subset = images(:, index);
subset_labels = labels(index);

% Store the subset in the data folder
csvwrite(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist' num2str(iter) '.csv'], subset);
csvwrite(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist' num2str(iter) 'labels.csv'], subset_labels);

cd(oldFolder);
end