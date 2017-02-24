function [] = preprocess_newsgroups( vocab_size )
% PREPROCESS_NEWSGROUPS Function to preprocess the 20 newsgroups dataset. 
% Creates a bag of words model with the top 'vocab_size' most frequent 
% words in the dataset.
%
%   INPUT:
%   vocab_size - number of words to keep in the vocabulary
%
%   OUTPUT:
%   Doesn't return anything. Saves the preprocessed data as a .mat file.
%

% Read the full vocabulary to get word indices
oldFolder = cd('C:\CMU\CMU-Spring-2016\DAP\Dataset\20newsgroups dataset');
fileID = fopen('vocabulary.txt');
full_vocab = textscan(fileID,'%s');
fclose(fileID);

% Read the train and test data and count the occurences of each word in the
% corpus
cd('20news-bydate\matlab')
train = load('train.data');
test = load('test.data');
data = [train; test];

vocab_count = zeros(size(full_vocab{1}, 1), 1);
for i = 1:size(data, 1)
    idx = data(i, 2);
    vocab_count(idx) = vocab_count(idx) + data(i, 3);
end

end









