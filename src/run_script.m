
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\iono-mod-norm.csv', 1, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\iono-mod-norm.csv', 2, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\iono-mod-norm.csv', 3, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\iono-mod-norm.csv', 4, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\iono-mod-norm.csv', 1, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\iono-mod-norm.csv', 2, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\iono-mod-norm.csv', 3, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\iono-mod-norm.csv', 4, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\iono-mod-norm.csv', 5, 2);
% 
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\CTG-mod-norm.csv', 1, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\CTG-mod-norm.csv', 2, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\CTG-mod-norm.csv', 3, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\CTG-mod-norm.csv', 4, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\CTG-mod-norm.csv', 1, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\CTG-mod-norm.csv', 2, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\CTG-mod-norm.csv', 3, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\CTG-mod-norm.csv', 4, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\CTG-mod-norm.csv', 5, 2);
% 
% 
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 1, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 2, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 3, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 4, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 1, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 2, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 3, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 4, 2);
main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 5, 2);
%
%
%
%
%
% D = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', 1, 1);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1_Dictionary_dp-max.mat', 'D');
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1 - dp-max algorithm');

% D = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', 2, 1);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1_Dictionary_dl-max.mat', 'D');
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1 - dl-max algorithm');

% D = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', 3, 1, 10);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1_Dictionary_dch-max.mat', 'D');
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1 - dch-max algorithm');

% D = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', 4, 1, 10);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1_Dictionary_dchperceptron-max.mat', 'D');
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1 - dchperceptron-max algorithm');

% D = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', 1, 2);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1_Dictionary_dp-mean.mat', 'D');
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1 - dp-mean algorithm');

% D = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', 2, 2);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1_Dictionary_dl-mean.mat', 'D');
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1 - dl-mean algorithm');

% D = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', 3, 2, 10);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1_Dictionary_dch-mean.mat', 'D');
% visualize_dict(D, 50, 28, 28, 'Dictionary atoms for mnist1 - dch-mean algorithm');

% D = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', 4, 2, 10);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1_Dictionary_dchperceptron-mean.mat', 'D');
% visualize_dict(D, 50, 28, 28, 'Dictionary atoms for mnist1 - dchperceptron-mean algorithm');

% D = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1.csv', 5, 2);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1_Dictionary_KSVD.mat', 'D');
% visualize_dict(D, 50, 28, 28, 'Dictionary atoms for mnist1 - KSVD algorithm');
%
%
%
%
%
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\20newsgroups1000_1.csv', 1, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\20newsgroups1000_1.csv', 2, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\20newsgroups1000_1.csv', 3, 1, 10);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\20newsgroups1000_1.csv', 4, 1, 10);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\20newsgroups1000_1.csv', 1, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\20newsgroups1000_1.csv', 2, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\20newsgroups1000_1.csv', 3, 2, 10);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\20newsgroups1000_1.csv', 4, 2, 10);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\20newsgroups1000_1.csv', 5, 2);


% % 
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist2.csv', 1, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist2.csv', 1, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist2.csv', 3, 1, 10);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist2.csv', 3, 2, 10);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist2.csv', 4, 1, 10);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist2.csv', 4, 2, 10);

% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist3.csv', 1, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist3.csv', 1, 2);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist4.csv', 1, 1);
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist4.csv', 1, 2);





%main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist2.csv', 5);
%main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist3.csv', 5);






