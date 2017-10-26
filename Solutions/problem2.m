image = imread('/Users/xuefanyong/Documents/MATLAB/skeleton_orig.tif');
[row,col] = size(image);
c_image = image;
d_image = image;
e_image = image;
f_image = image;
g_image = image;
h_image = image;
double_image = double(image);
lap_image = zeros(row,col);
mask = [0 1 0;1 -4 1;0 1 0];

%append_image = [zeros(row,1) double_image zeros(row,1)];
%append_image = [zeros(1,col+2); append_image ;zeros(1,col+2)];

%d
x_mask = [-1 -2 -1;0 0 0;1 2 1];
y_mask = [-1 0 1;-2 0 2;-1 0 1];

for r = 2:row-1
    for c = 2:col-1
        lap_image(r,c) = sum(sum(double_image(r-1:r+1,c-1:c+1).*mask));
        c_image(r,c) = image(r,c) - lap_image(r,c);
        d_image(r,c) = abs(sum(sum(double_image(r-1:r+1,c-1:c+1).*x_mask)))+abs(sum(sum(double_image(r-1:r+1,c-1:c+1).*y_mask)));
    end
end

%d_image

for r = 3:row-2
    for c = 3:col-2
        e_image(r,c) = mean(mean(double_image(r-2:r+2,c-2:c+2)));
    end
end

%f_image g_image h_image


f_image = e_image.*c_image;
g_image = f_image+image;
h_image = sqrt(double(g_image));


figure(1);
subplot(241);
imshow(image);
subplot(242);
imshow(lap_image);

subplot(245);
imshow(c_image);
subplot(246);
imshow(d_image);

subplot(243);
imshow(e_image);

subplot(244);
imshow(f_image);

subplot(247);
imshow(g_image);

subplot(248);
imshow(h_image);