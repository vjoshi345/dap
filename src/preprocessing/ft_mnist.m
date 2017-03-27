function [] = ft_mnist(data_path)
% FT_MNIST Function to store the 2-D Fourier Transform of mnist images
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
P_ft = zeros(d, n);
for i = 1:n
    image = reshape(P(:, i), [sz sz]);
    image_ft = fft2(image);
    P_ft(:, i) = reshape(image_ft, [d 1]);
    disp(['Iteration:' num2str(i) ' complete']);
end
dlmwrite(['C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\' file_name '-ft.csv'], P_ft);
end

