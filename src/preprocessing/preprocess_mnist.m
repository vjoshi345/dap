function [] = preprocess_mnist(dataset_size, iter)
% Script to preprocess and store the MNIST images in the data folder
% Stores a random subset of the images
rng(0);

% Read MNIST image data and get a random subset of it
oldFolder = cd('C:\CMU\CMU-Spring-2016\DAP\Dataset\MNIST database');
images = loadMNISTImages('train-images.idx3-ubyte');
[~, n] = size(images);
index = randi([1 n], 1, dataset_size);
subset = images(:, index);

% Store the subset in the data folder
csvwrite(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist' num2str(iter) '.csv'], subset);

cd(oldFolder);
end