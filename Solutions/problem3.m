%{
    Fourier Transform to be done
    Problem 3
    by Fanyong Xue 
    Student ID:515030910443
    Filtering in frequency domain
    
%}

%% Main Part
clear all
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/characters_test_pattern.tif');

ideal_low_plot_data(image);
ideal_high_plot_data(image);
butterworth_low_plot_data(image);
butterworth_high_plot_data(image);
gaussian_low_plot_data(image);
gaussian_high_plot_data(image);


%% Function Part

% Ideal
function ideal_high_plot_data(image)

    figure('name','Ideal High Pass');
    
    subplot(321);
    imshow(image);
    title('Original Image');
    
    subplot(322);
    imshow(IHPF(image,10),[]);
    title('Radius = 10');
    
    subplot(323);
    imshow(IHPF(image,30),[]);
    title('Radius = 30');
    
    subplot(324);
    imshow(IHPF(image,60),[]);
    title('Radius = 60');
    
    subplot(325);
    imshow(IHPF(image,160),[]);
    title('Radius = 160');
    
    subplot(326);
    imshow(IHPF(image,460),[]);
    title('Radius = 460');
    
end
function ideal_low_plot_data(image)

    figure('name','Ideal Low Pass');
    
    subplot(321);
    imshow(image);
    title('Original Image');
    
    subplot(322);
    imshow(ILPF(image,10),[]);
    title('Radius = 10');
    
    subplot(323);
    imshow(ILPF(image,30),[]);
    title('Radius = 30');
    
    subplot(324);
    imshow(ILPF(image,60),[]);
    title('Radius = 60');
    
    subplot(325);
    imshow(ILPF(image,160),[]);
    title('Radius = 160');
    
    subplot(326);
    imshow(ILPF(image,460),[]);
    title('Radius = 460');
    
end
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


% Butterworth
function butterworth_low_plot_data(image)
    figure('name','Butterworth Low Pass');
    
    subplot(321);
    imshow(image);
    title('Original Image');
    
    subplot(322);
    imshow(BLFP(image,2,10),[]);
    title('n=2,Radius = 10');
    
    subplot(323);
    imshow(BLFP(image,2,30),[]);
    title('n=2,Radius = 30');
    
    subplot(324);
    imshow(BLFP(image,2,60),[]);
    title('n=2,Radius = 60');
    
    subplot(325);
    imshow(BLFP(image,2,160),[]);
    title('n=2,Radius = 160');
    
    subplot(326);
    imshow(BLFP(image,2,460),[]);
    title('n=2,Radius = 460');
    
end
function butterworth_high_plot_data(image)
    figure('name','Butterworth High Pass');
    
    subplot(321);
    imshow(image);
    title('Original Image');
    
    subplot(322);
    imshow(BHFP(image,2,10),[]);
    title('n=2,Radius = 10');
    
    subplot(323);
    imshow(BHFP(image,2,30),[]);
    title('n=2,Radius = 30');
    
    subplot(324);
    imshow(BHFP(image,2,60),[]);
    title('n=2,Radius = 60');
    
    subplot(325);
    imshow(BHFP(image,2,160),[]);
    title('n=2,Radius = 160');
    
    subplot(326);
    imshow(BHFP(image,2,460),[]);
    title('n=2,Radius = 460');
    
end
function image_b = BLFP(image, n,radius)
    
    fliter = BLFP_fliter(image, n,radius);
    image_b = transfer(image,fliter);
    
end
function image_b = BHFP(image, n,radius)
    
    fliter = 1-BLFP_fliter(image, n,radius);
    image_b = transfer(image,fliter);
    
end
function fliter = BLFP_fliter(image, n,radius)
    [row,col]=size(image);
    fliter = zeros(2*row,2*col);
    
    for r = 1:2*row
        for c = 1:2*col
            fliter(r,c) = 1/(1+(sqrt((r-row)^2+(c-col)^2)/radius)^(2*n));
        end
    end
end

% Gaussian
function gaussian_low_plot_data(image)
    figure('name','Gaussian Low Pass');
    
    subplot(321);
    imshow(image);
    title('Original Image');
    
    subplot(322);
    imshow(GLFP(image,10),[]);
    title('Radius = 10');
    
    subplot(323);
    imshow(GLFP(image,30),[]);
    title('Radius = 30');
    
    subplot(324);
    imshow(GLFP(image,60),[]);
    title('Radius = 60');
    
    subplot(325);
    imshow(GLFP(image,160),[]);
    title('Radius = 160');
    
    subplot(326);
    imshow(GLFP(image,460),[]);
    title('Radius = 460');
end
function gaussian_high_plot_data(image)
    figure('name','Gaussian High Pass');
    
    subplot(321);
    imshow(image);
    title('Original Image');
    
    subplot(322);
    imshow(GHFP(image,10),[]);
    title('Radius = 10');
    
    subplot(323);
    imshow(GHFP(image,30),[]);
    title('Radius = 30');
    
    subplot(324);
    imshow(GHFP(image,60),[]);
    title('Radius = 60');
    
    subplot(325);
    imshow(GHFP(image,160),[]);
    title('Radius = 160');
    
    subplot(326);
    imshow(GHFP(image,460),[]);
    title('Radius = 460');
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
    
    for r = 1:2*row
        for c = 2:2*col
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
    image_t = abs(image_t);
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
%}