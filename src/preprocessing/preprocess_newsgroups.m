function [] = preprocess_newsgroups(vocab_size, ndocs, iter)
% PREPROCESS_NEWSGROUPS Function to preprocess the 20 newsgroups dataset. 
% Creates a bag of words model with the top 'vocab_size' most frequent 
% words in the dataset (stopwords are removed before computing frequency).
% Takes a random subset of the documents thus generated .
%
%   INPUT:
%   vocab_size - number of words to keep in the vocabulary
%   ndocs      - the size of random subset to be taken
%   iter       - an identifier for the size of subset chosen
%
%   OUTPUT:
%   Doesn't return anything. Saves the preprocessed data as a .csv file.
%   Saves the vocabulary as a .txt file.
%
%

rng(0);

% Read the full vocabulary to get word indices
oldFolder = cd('C:\CMU\CMU-Spring-2016\DAP\Dataset\20newsgroups dataset');
fileID = fopen('vocabulary.txt');
full_vocab = textscan(fileID,'%s');
full_vocab = full_vocab{1};
fclose(fileID);

% Remove stopwords
stopwords = {'a', 'about', 'above', 'above', 'across', 'after', ...
    'afterwards', 'again', 'against', 'all', 'almost', 'alone', 'along', ...
    'already', 'also','although','always','am','among', 'amongst', 'amoungst', ...
    'amount',  'an', 'and', 'another', 'any','anyhow','anyone','anything','anyway', ...
    'anywhere', 'are', 'around', 'as',  'at', 'back','be','became', 'because','become',...
    'becomes', 'becoming', 'been', 'before', 'beforehand', 'behind', 'being', 'below',...
    'beside', 'besides', 'between', 'beyond', 'bill', 'both', 'bottom','but', 'by',...
    'call', 'can', 'cannot', 'cant', 'co', 'con', 'could', 'couldnt', 'cry', 'de',...
    'describe', 'detail', 'do', 'done', 'down', 'due', 'during', 'each', 'eg', 'eight',...
    'either', 'eleven','else', 'elsewhere', 'empty', 'enough', 'etc', 'even', 'ever', ...
    'every', 'everyone', 'everything', 'everywhere', 'except', 'few', 'fifteen', 'fify',...
    'fill', 'find', 'fire', 'first', 'five', 'for', 'former', 'formerly', 'forty', 'found',...
    'four', 'from', 'front', 'full', 'further', 'get', 'give', 'go', 'had', 'has', 'hasnt',...
    'have', 'he', 'hence', 'her', 'here', 'hereafter', 'hereby', 'herein', 'hereupon', ...
    'hers', 'herself', 'him', 'himself', 'his', 'how', 'however', 'hundred', 'ie', 'if',...
    'in', 'inc', 'indeed', 'interest', 'into', 'is', 'it', 'its', 'itself', 'keep', 'last',...
    'latter', 'latterly', 'least', 'less', 'ltd', 'made', 'many', 'may', 'me', 'meanwhile',...
    'might', 'mill', 'mine', 'more', 'moreover', 'most', 'mostly', 'move', 'much', 'must',...
    'my', 'myself', 'name', 'namely', 'neither', 'never', 'nevertheless', 'next', 'nine',...
    'no', 'nobody', 'none', 'noone', 'nor', 'not', 'nothing', 'now', 'nowhere', 'of', 'off',...
    'often', 'on', 'once', 'one', 'only', 'onto', 'or', 'other', 'others', 'otherwise',...
    'our', 'ours', 'ourselves', 'out', 'over', 'own','part', 'per', 'perhaps', 'please',...
    'put', 'rather', 're', 'same', 'see', 'seem', 'seemed', 'seeming', 'seems', 'serious',...
    'several', 'she', 'should', 'show', 'side', 'since', 'sincere', 'six', 'sixty', 'so',...
    'some', 'somehow', 'someone', 'something', 'sometime', 'sometimes', 'somewhere', ...
    'still', 'such', 'system', 'take', 'ten', 'than', 'that', 'the', 'their', 'them',...
    'themselves', 'then', 'thence', 'there', 'thereafter', 'thereby', 'therefore', ...
    'therein', 'thereupon', 'these', 'they', 'thickv', 'thin', 'third', 'this', 'those',...
    'though', 'three', 'through', 'throughout', 'thru', 'thus', 'to', 'together', 'too',...
    'top', 'toward', 'towards', 'twelve', 'twenty', 'two', 'un', 'under', 'until', 'up',...
    'upon', 'us', 'very', 'via', 'was', 'we', 'well', 'were', 'what', 'whatever', 'when',...
    'whence', 'whenever', 'where', 'whereafter', 'whereas', 'whereby', 'wherein',...
    'whereupon', 'wherever', 'whether', 'which', 'while', 'whither', 'who', 'whoever',...
    'whole', 'whom', 'whose', 'why', 'will', 'with', 'within', 'without', 'would', 'yet',...
    'you', 'your', 'yours', 'yourself', 'yourselves', 'the'};
stopwords = stopwords';
rm_idx = ismember(full_vocab, stopwords);
all_indices = (1:size(full_vocab, 1))';
keep_idx = all_indices(~rm_idx);
full_vocab = full_vocab(keep_idx);

% Read the train and test data, remove stopwords, and count the occurences 
% of each word in the corpus
cd('20news-bydate\matlab')
train = load('train.data');
test = load('test.data');
data = [train; test];
data = full(sparse(data(:, 1), data(:, 2), data(:, 3)));
data = data(:, keep_idx);
vocab_count = sum(data);

% Sort the words by their frequency and keep only the top 'vocab_size' words
[~, sorted_idx] = sort(vocab_count, 'descend');
store_index = sorted_idx(1:vocab_size);
store_vocab = full_vocab(store_index);

cd('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data');
fid = fopen('newsgroups_vocab.txt', 'w');
[nrow, ~] = size(store_vocab);
formatSpec = '%s\n';
for row = 1:nrow
    fprintf(fid, formatSpec, store_vocab{row});
end
fclose(fid);

% Keep only the top 'vocab_size' words and then take a random subset of the
% data
data = data(:, store_index);
% index = randi([1 size(data, 1)], 1, ndocs);
[n, ~] = size(data);
index = datasample(1:n, ndocs, 'Replace', false);
subset = data(index, :);

% Normalize the dataset to have all feature values between 0 and 1
[~, d] = size(subset);
for i = 1:d
    if max(subset(:, i)) == min(subset(:, i))
        subset(:, i) = 1;
    else
        subset(:, i) = (subset(:, i) - min(subset(:, i)))/(max(subset(:, i)) - min(subset(:, i)));
    end
end
subset = subset';

csvwrite(['20newsgroups' num2str(vocab_size) '_' num2str(iter) '.csv'], subset);

cd(oldFolder);

end


















