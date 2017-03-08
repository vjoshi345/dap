function [] = visualize_dict(D, nimages, width, height, string)
% VISUALIZE_DICT Prints the dictionary atoms as images.
%   
%   INPUT:
%   D       - the dictionary to be visualized
%   nimages - number of atoms from the dictionary to be visualized (random
%             sample is taken)
%   width   - the width of the image to be plotted (aka #columns of the
%             matrix)
%   width   - the height of the image to be plotted (aka #rows of the
%             matrix)
%   string  - string to be put as the title
%
%   OUTPUT:
%   Doesn't output anything. The visualized image is saved.
%

rng(0);
assert(nimages <= size(D, 2), 'Cannot visualize more images than input!');
dict_to_visualize = datasample(D, nimages, 2, 'Replace', false);
assert(rem(nimages, 10) == 0, 'No. of images to print should be a multiple of 10!');
assert(width*height == size(D, 1), 'Height and width incorrect! Does not match size of input image!');

nrows = nimages/10;
figure();
for i = 1:nimages
    h = subplot(nrows, 10, i);
    p = get(h, 'pos');
    p(3) = p(3) + 0.015;
    p(4) = p(4) + 0.015;
    set(h, 'pos', p);
    curr = dict_to_visualize(:, i);
    curr = reshape(curr, [height, width]);
    imshow(curr);
end
suptitle(string);
end

