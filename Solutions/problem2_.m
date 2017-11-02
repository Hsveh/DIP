image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/skeleton_orig.tif');
[row,col] = size(image);
mask = [-1 -1 -1;-1 8 -1;-1 -1 -1];
mask = uint8(mask);
image_l = image;

for r = 2:row-1
    for c = 2:col-1
        image_l(r,c) = image(r,c)+sum(sum(image(r-1:r+1,c-1:c+1).*mask));
    end
end
%{
min = image_l(1,1);
max = image_l(1,1);

for r = 1:row
    for c = 1:col
        if min > image_l(r,c)
            min = image_l(r,c);
        end
        if max < image_l(r,c)
            max = image_l(r,c);
        end
    end
end

max = max - min;
image_l = image_l - min;
image_l = 150*(image_l/max);
image_l = image_l+100;
%}
figure();
subplot(242);
imshow(image_l);
subplot(243);
imshow(image_l+image);