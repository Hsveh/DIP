image = imread('/Users/xuefanyong/Documents/MATLAB/skeleton_orig.tif');
[row,col] = size(image);

%Fourier Transform
fourier_image = zeros(row,col);

%Transform Matrix
transform_r = zeros(row,1);
transform_c = zeros(col,1);
for r = 1:row
    transform_r(r) = exp(-j*2*pi*r/row);
end

for c = 1:col
    transform_c(c) = exp(-j*2*pi*c/col);
end

double_image = double(image);
for r = 1:row
    for c = 1:col
        %transform_m = transform_r.^r * (transform_c.^c)';
        %sum(sum(double_image.*transform_m));
        %fourier_image(r,c) = sum((double_image.*transform_m));
        %fourier_image(r,c) = fourier_image(r,c)/(row*col);
    end
end

%display all the images
figure(1);
subplot(2,4,1);
imshow(image);
title('Original Image');

subplot(2,4,2);
imshow(fourier_image);
title('Image after Fouier Transform');

subplot(2,4,3);
imshow(fft2(image));

