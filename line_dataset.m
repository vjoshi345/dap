% Script to create a dummmy dataset made up of points sampled from a set of
% lines

clear all;

% Setting the parameters
l = 10; % No. of lines
n = 100; % No. of points in each line
d = 10; % Dimension of each point

line_data = zeros(d, l*n);

rng(0);

for j = 1:l
  
    % Generate random numbers to define the line
    a = rand(1);
    b = rand(1);
    c = rand(1);
    
    % Generate random numbers for the rest of the points
    rest = rand(1, d-2);

    for i = 1:n
        r1 = rand(1);
        r2 = rand(1);

        temp = r1;

        r1 = r1*c/(a*r1 + b*r2);
        r2 = r2*c/(a*temp + b*r2);

        line_data(1, (j-1)*n + i) = r1;
        line_data(2, (j-1)*n + i) = r2;
        line_data(3:d, (j-1)*n + i) = rest;
    end
end

line_data = line_data';
dlmwrite([num2str(l) 'line-data-mod.csv'], line_data);
