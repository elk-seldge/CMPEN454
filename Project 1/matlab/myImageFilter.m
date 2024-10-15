function [img1] = myImageFilter(img0, h)
% input: original img, filter
% output: image after convolution

% pad: padarray()
% convolution
% for start row to end
%   for start col to end
%       for start row of filter to end
%           for start col of filter to end
[rows, cols] = size(img0); % img size
[fltrR, fltrC] = size(h);
pad_img0 = padarray(img0, [floor(fltrR / 2), floor(fltrC / 2)], "replicate");



[filterR, filterC] = size(h); % 3 * 3
convd = [(filterR - 1) / 2, (filterC - 1) / 2]; % 1 * 1

pad0 = padarray(img0, convd,"replicate"); % 17 * 17
img1 = zeros(rows, cols); % 16 *16
fltr = zeros(filterR, filterC); % 3 * 3

for i=1:filterR 
    for j=1:filterC
        fltr(i, j) = h(filterR - i + 1, filterC - j + 1); % flip the filter to prepare for conv
    end
end
for i = 1: rows
    for j = 1: cols
        temp = pad0(i:filterR + i - 1, j: filterC + j - 1) .* fltr; %cnvltn
        img1(i, j) = sum(temp(:)); % sum cnvltn val
    end
end

end

