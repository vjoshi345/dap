
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
% main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\wdbc-mod-norm.csv', 5, 2);
%
%
%
%
%
%
%
%
%

Y = csvread('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv');

% param.max_sparsity = 1; param.stopping_criterion = 1;
% out = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv', 1, param);
% idx = out.dict_idx;
% D = Y(:, idx);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1-ft_Dictionary_dp-max.mat', 'D');
% [d, n] = size(D);
% for i = 1:n
%     D(:, i) = reshape(ifft2(reshape(D(:, i), [28 28]), 'symmetric'), [d 1]);
% end
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1-ft - dp-max algorithm');

% param.max_sparsity = 2; param.stopping_criterion = 1;
% out = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv', 2, param);
% idx = out.dict_idx;
% D = Y(:, idx);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1-ft_Dictionary_dl-max.mat', 'D');
% [d, n] = size(D);
% for i = 1:n
%     D(:, i) = reshape(ifft2(reshape(D(:, i), [28 28]), 'symmetric'), [d 1]);
% end
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1-ft - dl-max algorithm');
% 
% param.max_sparsity = 10; param.stopping_criterion = 1;
% [out] = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv', 3, param);
% idx = out.dict_idx;
% D = Y(:, idx);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1-ft_Dictionary_dch-max.mat', 'D');
% [d, n] = size(D);
% for i = 1:n
%     D(:, i) = reshape(ifft2(reshape(D(:, i), [28 28]), 'symmetric'), [d 1]);
% end
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1-ft - dch-max algorithm');

% param.max_sparsity = 10; param.stopping_criterion = 1;
% out = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv', 4, param);
% idx = out.dict_idx;
% D = Y(:, idx);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1-ft_Dictionary_dchperceptron-max.mat', 'D');
% [d, n] = size(D);
% for i = 1:n
%     D(:, i) = reshape(ifft2(reshape(D(:, i), [28 28]), 'symmetric'), [d 1]);
% end
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1-ft - dchperceptron-max algorithm');

% param.max_sparsity = 1; param.stopping_criterion = 2;
% out = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv', 1, param);
% idx = out.dict_idx;
% D = Y(:, idx);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1-ft_Dictionary_dp-mean.mat', 'D');
% [d, n] = size(D);
% for i = 1:n
%     D(:, i) = reshape(ifft2(reshape(D(:, i), [28 28]), 'symmetric'), [d 1]);
% end
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1-ft - dp-mean algorithm');

% param.max_sparsity = 2; param.stopping_criterion = 2;
% out = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv', 2, param);
% idx = out.dict_idx;
% D = Y(:, idx);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1-ft_Dictionary_dl-mean.mat', 'D');
% [d, n] = size(D);
% for i = 1:n
%     D(:, i) = reshape(ifft2(reshape(D(:, i), [28 28]), 'symmetric'), [d 1]);
% end
% visualize_dict(D, 100, 28, 28, 'Dictionary atoms for mnist1-ft - dl-mean algorithm');

% param.max_sparsity = 10; param.stopping_criterion = 2;
% out = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv', 3, param);
% idx = out.dict_idx;
% D = Y(:, idx);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1-ft_Dictionary_dch-mean.mat', 'D');
% [d, n] = size(D);
% for i = 1:n
%     D(:, i) = reshape(ifft2(reshape(D(:, i), [28 28]), 'symmetric'), [d 1]);
% end
% visualize_dict(D, 20, 28, 28, 'Dictionary atoms for mnist1-ft - dch-mean algorithm');

% param.max_sparsity = 10; param.stopping_criterion = 2;
% out = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv', 4, param);
% idx = out.dict_idx;
% D = Y(:, idx);
% save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1-ft_Dictionary_dchperceptron-mean.mat', 'D');
% [d, n] = size(D);
% for i = 1:n
%     D(:, i) = reshape(ifft2(reshape(D(:, i), [28 28]), 'symmetric'), [d 1]);
% end
% visualize_dict(D, 20, 28, 28, 'Dictionary atoms for mnist1-ft - dchperceptron-mean algorithm');

param.dictionary_size = 50;
out = main('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist1-ft.csv', 5, param);
D = out.U;
save('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\output\mnist1-ft_Dictionary_KSVD.mat', 'D');
visualize_dict(D, 50, 28, 28, 'Dictionary atoms for mnist1-ft - KSVD algorithm');

%
%
%
%
%
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






