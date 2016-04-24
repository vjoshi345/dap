% Script to create a dummmy dataset made up of points sampled from a line
clear all;
l = 1; %No. of lines
n = 100; % No. of points
d = 10; % Dimension of each point

rng(0);

% Generate random numbers to define the line
a = rand(1);
b = rand(1);
c = rand(1);

rest = rand(1, d-2);
line_data = zeros(d, n);

for i = 1:n
    r1 = rand(1);
    r2 = rand(1);
    
    temp = r1;
    
    r1 = r1*c/(a*r1 + b*r2);
    r2 = r2*c/(a*temp + b*r2);
    
    line_data(1, i) = r1;
    line_data(2, i) = r2;
    line_data(3:d, i) = rest;
end
line_data = line_data';
dlmwrite('line-data-mod.csv', line_data);
