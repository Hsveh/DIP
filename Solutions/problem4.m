%{
    Problem 4
    by Fanyong Xue
	Student ID:515030910443
    Generating different types of noise and comparing different noise reduction methods
%}

%% Main Part
clear all
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/Fig0503.tif');
plot_data_noises(image);

circuit = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/Circuit.tif');
plot_data_noises(circuit);
plot_mean_filter(circuit);
plot_order_statistic_filter(circuit);

%% Function Part
% plot data by mean filters
function plot_mean_filter(image)
    figure('name','Mean Filters1');
    
    [row,col] = size(image);
    image = im2double(image);
    subplot(221);
    imshow(image);
    title('(a) X-ray image');
    
    n = uniform_noise(row,col,0,0.3);
    g = gaussian_noise(n,0,0.08);
    image_b = image+g;
    subplot(222);
    imshow(image_b);
    title('(b) Image corrupted by additive Gaussian noise');
    
    subplot(223);
    imshow(arithmetic_mean_filter(image_b,3,3));
    title('(c) Result of filtering with an arithmetic mean filter');
    
    subplot(224);
    imshow(geometric_mean_filter(image_b,3,3));
    title('(d) Result of filtering with a geometric mean filter');
    
    figure('name','Mean Filters2');
    
    subplot(221);
    image_a_ = impulse_noise(image,0.1,0,-1,0);
    imshow(image_a_);
    title('(a) Image corrupted by pepper noise');
    
    subplot(222);
    image_b_ = impulse_noise(image,0.1,0,1,0);
    imshow(image_b_);
    title('(b) Image corrupted by salt noise');
    
    
    subplot(223);
    imshow(contraharmonic_mean_filter(image_a_,1.5,3,3));
    title('(c) Result of filtering (a) with a contra-harmonic filter');
    
    subplot(224);
    imshow(contraharmonic_mean_filter(image_b_,-1.5,3,3));
    title('(c) Result of filtering (a) with a contra-harmonic filter');
    
    figure('name','Mean Filters3');
    
    subplot(121);
    imshow(contraharmonic_mean_filter(image_a_,-1.5,3,3));
    title('(c) Result of filtering (a) with a contra-harmonic filter(Q=-1.5)');
    
    subplot(122);
    imshow(contraharmonic_mean_filter(image_b_,1.5,3,3));
    title('(c) Result of filtering (a) with a contra-harmonic filter(Q=1.5)');
end
% plot data by order statistic filters
function plot_order_statistic_filter(image)
    figure('name','Order-Statistic Filters1');
    image = im2double(image);
    image_a = impulse_noise(image,0.1,0.1,-1,1);
    
    subplot(221);
    imshow(image_a);
    title('Image corrupted by salt-and-pepper noise');
    
    subplot(222);
    image_b = median_filter(image_a,3,3);
    imshow(image_b);
    title('(b) Result of one pass with a median filter');
    
    subplot(223);
    image_c = median_filter(image_b,3,3);
    imshow(image_c);
    title('(c) Result of processing (b) with the same filter');
    
    subplot(224);
    imshow(median_filter(image_c,3,3));
    title('(d) Result of processing (c) with the same filter');
    
    figure('name','Order-Statistic Filters2');
    
    image_a__ = impulse_noise(image,0.1,0,-1,0);
    subplot(121);
    imshow(max_filter(image_a__,3,3));
    title('(a) Result of filtering pepper noise with a max filter');
    
    image_b__ = impulse_noise(image,0.1,0,1,0);
    subplot(122);
    imshow(min_filter(image_b__,3,3));
    title('(a) Result of filtering salt noise with a min filter');
    
    figure('name','Order-Statistic Filters3');
    [row,col] = size(image);
    n = uniform_noise(row,col,0,0.3);
    image_a_ = image + n;
    subplot(321);
    imshow(image_a_);
    title('(a) Image corrupted by additive uniform noise');
    
    image_b_ = impulse_noise(image_a_,0.1,0.1,-1,1);
    subplot(322);
    imshow(image_b_);
    title('(b) Image corrupted by additive salt-and-pepper noise');
    
    subplot(323);
    imshow(arithmetic_mean_filter(image_b_,5,5));
    title('(c) Image (b) filtered with an arithmetic mean filter');
    
    subplot(324);
    imshow(geometric_mean_filter(image_b_,5,5));
    title('(c) Image (b) filtered with a geometric mean filter');
    
    subplot(325);
    imshow(median_filter(image_b_,5,5));
    title('(c) Image (b) filtered with a median filter');
    
    subplot(326);
    imshow(alpha_trimmed_mean_filter(image_b_,5,5,5));
    title('(c) Image (b) filtered with an alpha-trimmed mean filter');
end
% arithmetic mean filter with m*n mask
function image_ = arithmetic_mean_filter(image,m,n)
    [row,col] = size(image);
    image_ = image;
    %m and n should be odd numbers
    m_ = floor(m/2);
    n_ = floor(n/2);
    for r = ceil(m/2):row-m_
        for c = ceil(n/2):col-n_
            image_(r,c) = mean(mean(image(r-m_:r+m_,c-n_:c+n_)));
        end
    end
end
% geometric mean filter with m*n mask
function image_ = geometric_mean_filter(image,m,n)
    [row,col] = size(image);
    image_ = image;
    m_ = floor(m/2);
    n_ = floor(n/2);
    for r = ceil(m/2):row-m_
        for c = ceil(n/2):col-n_
            image_(r,c) = nthroot(prod(prod(image(r-m_:r+m_,c-n_:c+n_))),m*n);
        end
    end
end
% harmonic mean filter with m*n mask
function image_ = harmonic_mean_filter(image,m,n)
    [row,col] = size(image);
    image_ = image;
    m_ = floor(m/2);
    n_ = floor(n/2);
    for r = ceil(m/2):row-m_
        for c = ceil(n/2):col-n_
            image_(r,c) = (m*n)/(sum(sum(1./image(r-m_:r+m_,c-n_:c+n_))));
        end
    end
end
% contraharmonic mean filter with m*n mask and its oder is q
function image_ = contraharmonic_mean_filter(image,q,m,n)
    [row,col] = size(image);
    image_ = image;
    m_ = floor(m/2);
    n_ = floor(n/2);
    for r = ceil(m/2):row-m_
        for c = ceil(n/2):col-n_
            image_(r,c) = (sum(sum(image(r-m_:r+m_,c-n_:c+n_).^(q+1))))/(sum(sum(image(r-m_:r+m_,c-n_:c+n_).^q)));
        end
    end
    image_ = real(image_);
end
% median filter with m*n mask
function image_ = median_filter(image,m,n)
    [row,col] = size(image);
    image_ = image;
    m_ = floor(m/2);
    n_ = floor(n/2);
    for r = ceil(m/2):row-m_
        for c = ceil(n/2):col-n_
            image_(r,c) = median(median(image(r-m_:r+m_,c-n_:c+n_)));
        end
    end
end
% max filter with m*n mask
function image_ = max_filter(image,m,n)
    [row,col] = size(image);
    image_ = image;
    m_ = floor(m/2);
    n_ = floor(n/2);
    for r = ceil(m/2):row-m_
        for c = ceil(n/2):col-n_
            temp = image(r-m_:r+m_,c-n_:c+n_);
            image_(r,c) = max(temp(:));
        end
    end
end
% min filter with m*n mask
function image_ = min_filter(image,m,n)
    [row,col] = size(image);
    image_ = image;
    m_ = floor(m/2);
    n_ = floor(n/2);
    for r = ceil(m/2):row-m_
        for c = ceil(n/2):col-n_
            temp = image(r-m_:r+m_,c-n_:c+n_);
            image_(r,c) = min(temp(:));
        end
    end
end
% midpoint filter with m*n mask
function image_ = midpoint_filter(image,m,n)
    [row,col] = size(image);
    image_ = image;
    m_ = floor(m/2);
    n_ = floor(n/2);
    for r = ceil(m/2):row-m_
        for c = ceil(n/2):col-n_
            temp = image(r-m_:r+m_,c-n_:c+n_);
            image_(r,c) = (max(temp(:))+min(temp(:)))/2;
        end
    end
end
% alpha trimmed mean filter with m*n mask(deleta d pixels)
function image_ = alpha_trimmed_mean_filter(image,d,m,n)
    [row,col] = size(image);
    image_ = image;
    m_ = floor(m/2);
    n_ = floor(n/2);
    for r = ceil(m/2):row-m_
        for c = ceil(n/2):col-n_
            temp = image(r-m_:r+m_,c-n_:c+n_);
            temp_ = sort(temp(:));
            image_(r,c) = mean(temp_(floor(d/2):m*n-floor(d/2)));
        end
    end
end

% adding noises to image and plot them
function plot_data_noises(image)

    [row,col] = size(image);
    image = im2double(image);
    figure('name','Uniform Noise');
    n = uniform_noise(row,col,0,0.3);
    image_=image+n;
    subplot(121);
    imshow(image_,[]);
    subplot(122);
    bar(get_histogram(image_));

    figure('name','Gaussian Noise');
    g = gaussian_noise(n,0,0.08);
    image_ = image+g;
    subplot(121);
    imshow(image_,[]);
    subplot(122);
    bar(get_histogram(image_));

    figure('name','Rayleigh Noise');
    r = rayleigh_noise(n,-0.2,0.03);
    image_ = image+r;
    subplot(121);
    imshow(image_,[]);
    subplot(122);
    bar(get_histogram(image_));

    
    figure('name','Exponential Noise');
    e = exponential_noise(n,25);
    image_ = image+e;
    subplot(121);
    imshow(image_,[]);
    subplot(122);
    bar(get_histogram(image_));

    figure('name','Gamma Noise');
    ga = gamma_noise(n,25,3);
    image_ = image+ga;
    subplot(121);
    imshow(image_,[]);
    subplot(122);
    bar(get_histogram(image_));

    figure('name','Impulse Noise');
    image_ = impulse_noise(image,0.1,0.1,0.2,-0.2);
    subplot(121);
    imshow(image_,[]);
    subplot(122);
    bar(get_histogram(image_));
end
% row*col uniform noise from low~high
function n = uniform_noise(row,col,low,high)
    n = low + (high-low)*rand([row col]);
end
% normalizing uniform noise to 0~1
function n = normalized(uniform_noise)
    max_ = max(uniform_noise(:));
    min_ = min(uniform_noise(:));
    n = double(uniform_noise-min_);
    n = n/double(max_-min_);
end
% rayleigh noise
function n = rayleigh_noise(uniform_noise,a,b)
    uniform_noise_ = normalized(uniform_noise);
    n = a + sqrt(-b*log(1-uniform_noise_));
end
% exponential noise
function n = exponential_noise(uniform_noise,a)
    uniform_noise_ = normalized(uniform_noise);
    n = -1/a*log(1-uniform_noise_);
end
% impulse noise
function image_ = impulse_noise(image,pa,pb,a,b)
    [row,col]=size(image);
    uniform_noise_ = rand([row col]);
    image_ = image;
    for r = 1:row
        for c = 1:col
            if uniform_noise_(r,c)<pa
                image_(r,c) = image(r,c)+a;
            elseif uniform_noise_(r,c)>(1-pb)
                image_(r,c) = image(r,c)+b;
            end
        end
    end
end
% gamma noise
function n = gamma_noise(uniform_noise_,a,b)
    uniform_noise_ = normalized(uniform_noise_);
    n = -1/a*log(1-uniform_noise_);
    [row,col] = size(uniform_noise_);
    for i = 2:b
        uniform_noise__=uniform_noise(row,col,0,1);
        n = n-1/a*log(1-uniform_noise__);
    end
    
end
% gaussian noise
function n = gaussian_noise(uniform_noise,ex,sigma)
    [row,col] = size(uniform_noise);
    n_ = normalized(uniform_noise);
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
% plot image' histogram (but add other 150 to display the negative value)
function histogram = get_histogram(image)
    [row,col]=size(image);
    histogram = zeros(500,1);%-150~350
    for r = 1:row
        for c = 1:col
            gray = int16(image(r,c) /0.004);
            if gray >349
                gray = 349;
            end
            if gray <-150
                gray = -150;
            end
            histogram(gray+1+150)=histogram(gray+1+150)+1;
        end
    end
end