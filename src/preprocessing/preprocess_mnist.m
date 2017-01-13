% Script to preprocess and store the MNIST images in the data folder
% Stores a random subset of the images
rng(0);

% Read MNIST image data and get a random subset of it
oldFolder = cd('C:\CMU\CMU-Spring-2016\DAP\Dataset\MNIST database');
images = loadMNISTImages('train-images.idx3-ubyte');
[~, n] = size(images);
index = randi([1 n], 1, 1000);
subset = images(:, index);

% Store the subset in the data folder
csvwrite('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', subset);

cd(oldFolder);