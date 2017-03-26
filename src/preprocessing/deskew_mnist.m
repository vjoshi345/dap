function [] = deskew_mnist(data_path)
% DESKEW_MNIST Preprocess mnist data to straighten the slanting digits.
% Note: Uses the horizon function downloaded from: https://www.mathworks.com/matlabcentral/fileexchange/40239-straighten-image
%
%   INPUT:
%   data_path - full path to the mnist data 
%
%   OUTPUT:
%   No variables returned. Saves the modified dataset at the default
%   location.
%
P = csvread(data_path);
[~, file_name, ~] = fileparts(data_path);
[d, n] = size(P);
sz = sqrt(d);
P_deskewed = zeros(d, n);
for i = 1:n
    image = reshape(P(:, i), [sz sz]);
    angle = horizon(image, 0.1, 'fft');
    P_deskewed(:, i) = reshape(imrotate(image, -angle, 'bicubic', 'crop'), [d 1]);
    disp(['Iter:' num2str(i) ', angle:' num2str(angle)]);
end
dlmwrite(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\' file_name '-deskewed1.csv'], P_deskewed);
end

