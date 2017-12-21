% image append and bugs to be solve
% Problem 2
% by Fanyong Xue
% Student ID:515030910443
% Combining spatial enhancement methods

%% Main Part
clear all
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/skeleton_orig.tif');
[row,col] = size(image);
mask = [-1 -1 -1;-1 8 -1;-1 -1 -1];
mask = double(mask);
b_image = laplace_transformations(image,mask);
figure();
imshow(b_image,[])
title('b');
c_image = b_image+im2double(image);
figure();
imshow(c_image,[])
title('c');
d_image = sobel_gradient(image);
figure();
imshow(d_image)
title('d');
e_image = smooth(d_image);
figure();
imshow(e_image)
title('e');
f_image = im2double(e_image).*c_image;
figure();
imshow(f_image,[])
title('f');
g_image = abs(f_image)+im2double(image);
figure();
imshow(g_image,[])
title('g');
h_image = sqrt(g_image);
figure();
imshow(h_image,[])
title('h');
plot_data(image,b_image,c_image,d_image,e_image,f_image,g_image,h_image);

%% Function Part

% Laplace Transfromation for image using mask
% Input:
%   image:image you want to perform
%   mask:Laplace mask you want to use
% Output:
%   image_l:image after laplace transformation

function image_l = laplace_transformations(image,mask)
    [row,col] = size(image);
    mask = double(mask);
    %append image
    image_l = im2double(image);
    image = [zeros(row,2) image zeros(row,2)];
    image = [zeros(2,col+4);image;zeros(2,col+4)];
    image_append = im2double(image);
    
    for r = 1:row
        for c = 1:col
            image_l(r,c) = sum(sum(image_append(r:r+2,c:c+2).*mask));
        end
    end
end
%{
    sobel gradient for image
%}
function image_s = sobel_gradient(image)
    [row,col] = size(image);
    x_mask = [-1 -2 -1;0 0 0;1 2 1];
    y_mask = [-1 0 1;-2 0 2;-1 0 1];
    image_s = image;
    image = double(image);
    
    for r = 2:row-1
        for c = 2:col-1
            image_s(r,c) = abs(sum(sum(image(r-1:r+1,c-1:c+1).*x_mask)))+abs(sum(sum(image(r-1:r+1,c-1:c+1).*y_mask)));
        end
    end
end
%{
    smooth image using 5*5 mean filter
%}
function image_s = smooth(image)
    [row,col] = size(image);
    image_s = image;
    for r = 3:row-2
        for c = 3:col-2
            image_s(r,c) = mean(mean(image(r-2:r+2,c-2:c+2)));
        end
    end
end
%{
    plot data
%}
function plot_data(a,b,c,d,e,f,g,h)
    figure();
    
    subplot(241);
    imshow(a);
    title('(a) Oringinal Image');
    
    subplot(242);
    imshow(b,[]);
    title('(b) Laplacian of (a)');
    
    subplot(245);
    imshow(c,[]);
    title('(c) Sharpened image');
    
    subplot(246);
    imshow(d);
    title('(d) Sobel gradient');
    
    subplot(243);
    imshow(e);
    title('(e) Smoothed sobel image');
    
    subplot(244);
    imshow(f,[]);
    title('(f) Product of (c) and (e)');
    
    subplot(247);
    imshow(g,[]);
    title('(g) Sharpened image');
    
    subplot(248);
    imshow(h,[]);
    title('(h) Final result');
end