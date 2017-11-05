%{
%}
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/Fig0503.tif');
plot_data_noises(image);

circuit = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/Circuit.tif');
plot_data_noises(circuit);


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
end
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
    im = impulse_noise(n,0.15);
    image_ = image+im;
    subplot(121);
    imshow(image_,[]);
    subplot(122);
    bar(get_histogram(image_));
end
function n = uniform_noise(row,col,low,high)
    n = low + (high-low)*rand([row col]);
end
function n = normalized(uniform_noise)
    max_ = max(uniform_noise(:));
    min_ = min(uniform_noise(:));
    n = double(uniform_noise-min_);
    n = n/double(max_-min_);
end
function n = rayleigh_noise(uniform_noise,a,b)
    uniform_noise_ = normalized(uniform_noise);
    n = a + sqrt(-b*log(1-uniform_noise_));
end
function n = exponential_noise(uniform_noise,a)
    uniform_noise_ = normalized(uniform_noise);
    n = -1/a*log(1-uniform_noise_);
end
function n = impulse_noise(uniform_noise,threshold)
    [row,col]=size(uniform_noise);
    n = double(zeros(row,col));
    for r = 1:row
        for c = 1:col
            if uniform_noise(r,c)>threshold
                n(r,c) = 1;
            else
                n(r,c) = 0;
            end
        end
    end
end
function n = gamma_noise(uniform_noise_,a,b)
    uniform_noise_ = normalized(uniform_noise_);
    n = -1/a*log(1-uniform_noise_);
    [row,col] = size(uniform_noise_);
    for i = 2:b
        uniform_noise__=uniform_noise(row,col,0,1);
        n = n-1/a*log(1-uniform_noise__);
    end
    
end
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