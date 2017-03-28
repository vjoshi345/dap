load('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\to_publish\Visualizations\Dictionary\Deskewed\mnist1-deskewed_Dictionary_dch-max.mat');
Y = csvread('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-deskewed.csv');
param.niter = 10;
param.epsilon = 4.66093;
blah = setdiff(Y', D', 'rows')';
visualize_reconstruction(D, blah(:, 10), 1, param);
visualize_reconstruction(D, blah(:, [50, 73]), 1, param);
visualize_reconstruction(D, blah(:, [100, 158, 366, 513]), 1, param);