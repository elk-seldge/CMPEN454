function [rhos, thetas] = myHoughLines(H, nLines)
%Your implemention here
%nLines = 30;
% first, NMS all the neighbours;
% Then use builtin sort() to sort
% rhos = zeros(nLines, 1);
% thetas = zeros(nLines, 1);
% while (index < nLine)
% find corresponding index
% rhos(index, 1); thetas(index, 1); index++;

rhos = zeros(nLines, 1);
thetas = zeros(nLines, 1);
nmsResult = H;
[rows, cols] = size(nmsResult);
HPadded = padarray(nmsResult, [1, 1]);

for i = 2: rows - 1 % NMS part
    for j = 2: cols - 1 % any() check if there is any non-zero value in matrix
        if any(find(HPadded(i - 1: i + 1, j - 1: j + 1) > HPadded(i, j))) % compare with neighbours
            nmsResult(i, j) = 0;
        end
    end
end
nmsResult = nmsResult(2:rows-1,2:cols-1);
sorted = sort(nmsResult,'descend');
sorted = sorted(1:nLines, 1);
num = 0;
while num < nLines
    num = num + 1;
    [r, c] = find(nmsResult == sorted(num));
    nmsResult(r(1), c(1))=0;
    rhos(num,1) = r(1);
    thetas(num,1) = c(1);    
end
end
        