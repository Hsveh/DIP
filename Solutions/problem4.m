%{
%}
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/Fig0503.tif');

[row,col] = size(image);

figure();
n = uniform_noise(row,col,0,50);
k=image+n;
subplot(221);
imshow(k);
subplot(222);
imhist(k);
g=impulse_noise(n,25);
kk = image + g;
subplot(223);
imshow(kk);
subplot(224);
imhist(kk);
function n = uniform_noise(row,col,low,high)
    n = uint8(low + (high-low)*rand([row col]));
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
    n = uint8(n);
end
function n = exponential_noise(uniform_noise,a)
    uniform_noise_ = normalized(uniform_noise);
    n = -1/a*log(1-uniform_noise_);
end
function n = impulse_noise(uniform_noise,threshold)
    [row,col]=size(uniform_noise);
    n = uint8(zeros(row,col));
    for r = 1:row
        for c = 1:col
            if uniform_noise(r,c)>threshold
                n(r,c) = 255;
            else
                n(r,c) = 0;
            end
        end
    end
end
function n = gamma_noise(uniform_noise,a,b)
    uniform_noise_ = normalized(uniform_noise);
    
end

%{
function n = gaussian_noise(uniform_noise,ex,sigma)
    distribution = zeros(256,1);
    distribution(1) = gaussian(1,ex,sigma);
    
    for z = 2:256
        distribution(z) = distribution(z-1) + gaussian(z,ex,sigma);
    end
    
    max_ = max(uniform_noise(:));
    min_ = min(uniform_noise(:));
    
    mapping = zeros(256,1);
    
    for z = 1:256
        mapping(z) = min_+round((max_-min_)*distribution(z));
    end
    
    mapping_re = zeros(max_-min_+1,1);
    mapping_re(1) = 0;
    for z = 2:256
        if ~(mapping(z-1)+1<mapping(z))
            continue;
        end
        for z_ = mapping(z-1)+1 : mapping(z)
            mapping_re(z_+1-min_) = z-1;
        end
    end
    
    [row,col] = size(uniform_noise);
    n = uint8(zeros(row,col));
    for r = 1:row
        for c = 1:col
            n(r,c) = mapping_re(uniform_noise(r,c)+1-min_);
        end
    end
end
%}

function n = gaussian_noise(uniform_noise,ex,sigma)
    n_ = uniform_noise;
    [row,col] = size(uniform_noise);
    
    max_ = max(uniform_noise(:));
    min_ = min(uniform_noise(:));
    n_ = n_-min_;
    n_ = double(n_)/double(max_-min_);
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
    n=uint8(n);
end
