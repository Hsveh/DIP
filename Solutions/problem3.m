image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/characters_test_pattern.tif');
[row,col] = size(image);

%Fourier Transform
%fourier_image = zeros(row,col);

%Transform Matrix
%transform_r = zeros(row,1);
%transform_c = zeros(col,1);
%for r = 1:row
%    transform_r(r) = exp(-j*2*pi*r/row);
%end

%for c = 1:col
%    transform_c(c) = exp(-j*2*pi*c/col);
%end

%double_image = double(image);
%for r = 1:row
%    for c = 1:col
        %transform_m = transform_r.^r * (transform_c.^c)';
        %sum(sum(double_image.*transform_m));
        %fourier_image(r,c) = sum((double_image.*transform_m));
        %fourier_image(r,c) = fourier_image(r,c)/(row*col);
%    end
%end

%display all the images
figure(1);
subplot(2,4,1);
imshow(image);
title('Original Image');

subplot(2,4,2);


imshow(ILPF(image,60),[]);
title('Image after Fouier Transform');


function image_i = ILPF(image,radius)
    
    fliter = ILPF_fliter(image,radius);
    image_i = transfer(image,fliter);
    %{
    image_f = fft2(image,2*row,2*col);
    image_f = fftshift(image_f);
    %image_f = log(1+abs(image_f));
    
    %%%%%%%
    image_i = image_f.*fliter;
    
    image_i = ifftshift(image_i);
    image_i = ifft2(image_i);
    
    %image_i = abs(image_i);
    image_i = image_i(1:row,1:col);
    image_i = real(image_i);
    image_f = log(1+abs(image_f));
    %}
end

function image_i = IHPF(image,radius)
    fliter = 1-ILPF_fliter(image,radius);
    image_i = transfer(image,fliter);
end

function fliter = ILPF_fliter(image,radius)
    [row,col]=size(image);
    fliter = zeros(2*row,2*col);
    
    for r = 1:2*row
        for c = 1:2*col
            if sqrt((r-row)^2+(c-col)^2) <= radius
                fliter(r,c) = 1;
            end
        end
    end
end

function image_b = BLFP(image, n,radius)
    
    fliter = BLFP_fliter(image, n,radius);
    image_b = transfer(image,fliter);
    
end

function image_b = BGFP(image, n,radius)
    
    fliter = 1-BLFP_fliter(image, n,radius);
    image_b = transfer(image,fliter);
    
end

function fliter = BLFP_fliter(image, n,radius)
    [row,col]=size(image);
    fliter = double(zeros(2*row,2*col));
    
    for r = 1:row
        for c = 2:col
            fliter(r,c) = 1/(1+(sqrt((r-row)^2+(c-col)^2)/radius)^(2*n));
        end
    end
end

function image_g = GLFP(image,radius)
    fliter = GLFP_fliter(image,radius);
    image_g = transfer(image,fliter);
end

function image_g = GHFP(image,radius)
    fliter = 1-GLFP_fliter(image,radius);
    image_g = transfer(image,fliter);
end

function fliter = GLFP_fliter(image,radius)
    [row,col]=size(image);
    fliter = double(zeros(2*row,2*col));
    
    for r = 1:row
        for c = 2:col
            fliter(r,c) = exp(-1*(((r-row)^2+(c-col)^2)/(2*radius^2)));
        end
    end
end

function image_t=transfer(image,fliter)
    [row,col]=size(image);
    image_f = fft2(image,2*row,2*col);
    image_f = fftshift(image_f);
    image_t = image_f.*fliter;
    image_t = ifftshift(image_t);
    image_t = ifft2(image_t);
    image_t = image_t(1:row,1:col);
    image_t = real(image_t);
end
%{

function [ mfft2 ] = JCGuoFFT2( data )
    h = size(data, 1);
    w = size(data, 2);
    mfft2 = data;

    if power(2, log2(h)) ~= h || power(2, log2(w)) ~= w
        disp('JCGuoFFT2 exit: h and w must be the power of 2!')
    else
        for i = 1 : h
            mfft2(i, :) = IterativeFFT(mfft2(i, :));
        end

        for j = 1 : w
            mfft2(:, j) = IterativeFFT(mfft2(:, j));
        end
    end
end

function image_s = shift_image(image)

    [row,col]=size(image);
    image_s = image;
    
    for r = 1:row
        for c = 1:col
            image_s(r,c) = image(r,c)*(-1)^(r+c);
        end
    end
end


function image_f = DFT(image,rows,cols)
    [row,col]=size(image);
    
    %pad image to rows*cols
    image = [image zeros(row,cols-col)];
    image = [image;zeros(rows-row,cols)];
    image = double(image);
    for i = 1:rows
        k = cols/2;
        M = round(log2(cols));
        for j = 1:cols-2
            if j<k
                t = image(i,k);
                image(i,k) = image(i,j);
                image(i,j) = t;
            end
            l = cols/2;
            while l<=k
                k = k-1;
                l = l/2;
            end
            k = k+1;
        end
        for m = 1:M
            la = 2^m;
            lb = la/2;
            for l = 1:lb
                r = (l-1)*2^(M-m);
                n = l-1;
                while n<rows-1
                    lc = n+lb;
                    t = image(lc,j)*exp(2*pi*r/rows);
                    image(i,lc) = image(i,n) - t;
                    image(i,n) = image(i,n) + t;
                    n = n+la;
                end
            end
        end
    end
    image_f = image;
end
%}

function image_f=DFT(image,rows,cols)

end

function v = DFT_1(V)
    n = length(V);
    fft_m = BitReverseCopy(V);
    
    for r = 1:log2(n)
        m = power(2,r);
        wm = exp(- 2 * pi * i / m);
        
        for k = 0 : m : n - 1
            w = 1;
            for j = 0 : m / 2 - 1
                t = w * fft_m(k + j + m / 2 + 1);
                u = fft_m(k + j + 1);
                fft_m(k + j + 1) = u + t;
                fft_m(k + j + m / 2 + 1) = u - t;
                w = w * wm;
            end
        end
    end
end