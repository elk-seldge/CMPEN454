gaus = [1, 2, 1;
        2, 4, 2;
        1, 2, 3;];
gsize = size(gaus);
img_test = [1, 2, 3, 4;
            5, 7, 11, 2;
            3, 4, 5, 56;];
[rows, cols] = size(img_test);
pad_img_test = padarray(img_test, [floor(gsize(1)/2), floor(gsize(2)/2)],"replicate");
flipud(im2col(gaus, [gsize(1), gsize(2)]))
fltr = fliplr(flipud(gaus));
img1 = zeros(rows, cols);
for i=1:rows
    for j=1:cols
        temp = pad_img_test(i:i+gsize(1)-1, j:j+gsize(2)-1) .* fltr;
        img1(i, j) = sum(temp(:));
    end
end
img1
xcorr(fltr, img_test)
xcorr(img_test, fltr)
