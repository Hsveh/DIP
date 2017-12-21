%{
    Problem 5
    by Fanyong Xue
	Student ID:515030910443
    Image restoration
%}
%% Main Part

clear all
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/book_cover.jpg');
image = im2double(image);

% motion blur
[row,col] = size(image);
image_f = fft2(image,row,col);
image_f2 = fftshift(image_f);
fliter = blurring_degradation(row,col,0.1,0.1,1);
image_t = image_f2.*fliter;
G = image_t;
image_t = ifftshift(image_t);
image_t = ifft2(image_t);
image_t = abs(image_t);
image_t = im2uint8(image_t);

% add noise
noise=gaussian_noise(row,col,0,25);
noise=int16(noise);
image_n = noise+int16(image_t);
image_n = uint8(image_n);
noise=uint8(noise);

% diaplsy
figure();
subplot(221);
imshow(image);
title('Original Image');
subplot(222);
imshow(image_t);
title('Blurred Image');
subplot(223);
imshow(noise);
title('Noise');
subplot(224);
imshow(image_n);
title('Blurred Image with Noise');

% restoration
% no noise
image_i = inverse_filter(G,fliter);
figure();
imshow(image_i);
title('Inverse Filter');

% with noise 650
N = fft2(im2double(noise));
N = fftshift(N);
G_ = G+N;
image_i2 = inverse_filter(G_,fliter);
image_w = wiener_filter(G_,fliter,0.03);
% display
figure();
subplot(121);
imshow(image_i2,[]);
title('Inverse Filter(650)');
subplot(122);
imshow(image_w);
title('Wiener Filter(650)');

%with noise 324
noise=gaussian_noise(row,col,0,18);
noise=uint8(noise);
N = fft2(im2double(noise));
N = fftshift(N);
G_ = G+N;

image_i2 = inverse_filter(G_,fliter);
image_w = wiener_filter(G_,fliter,0.01);

% display
figure();
subplot(121);
imshow(image_i2,[]);
title('Inverse Filter(324)');
subplot(122);
imshow(image_w);
title('Wiener Filter(324)');
% 
%% Function Part
% blurring degradation
function h = blurring_degradation(row,col,a,b,T)
    
    h = double(zeros(row,col));
    for u  = 1:row
        for v = 1:col
            if u +v == col
                h(u,v)=1;
                continue;
            end
            h(u,v) = T/(pi*((u-row/2)*a+(v-col/2)*b))*sin(pi*((u-row/2)*a+(v-col/2)*b))*exp(-1i*pi*((u-row/2)*a+(v-col/2)*b));
        end
    end
end

% gaussian noise
function n = gaussian_noise(row,col,ex,sigma)
    n_ = rand([row col]);
    n = n_;
    for r = 1:row
        if mod(col,2)==0
            for c = 1:2:col
                n(r,c) = sqrt(-2*log(n_(r,c)))*cos(2*pi*n_(r,c+1));
                n(r,c+1) = sqrt(-2*log(n_(r,c)))*sin(2*pi*n_(r,c+1));
            end
        else
            for c = 1:2:col-1
                n(r,c) = sqrt(-2*log(n_(r,c)))*cos(2*pi*n_(r,c+1));
                n(r,c+1) = sqrt(-2*log(n_(r,c)))*sin(2*pi*n_(r,c+1));
            end
            n(r,col) = sqrt(-2*log(n_(r,c)))*cos(2*pi*n_(r,1));
        end
    end
    n = n*double(sigma);
    n = n+ex;
end
% inverse filter
function image_ = inverse_filter(g,h)
    image_ = g./h;
    image_ = ifftshift(image_);
    image_ = ifft2(image_);
    image_ = abs(image_);
end
% wiener filter
function image_ = wiener_filter(g,h,k)
    image_ = ((conj(h).*h)./((conj(h).*h)+k)).*g./h;
    image_ = ifftshift(image_);
    image_ = ifft2(image_);
    image_ = abs(image_);
end