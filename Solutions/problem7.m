%{
%}
clear all
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/lenna.tif');
image = im2double(image);
zonal_mask = [1 1 1 1 1 ;1 1 1 1 0;1 1 1 0 0;1 1 0 0 0;1 0 0 0 0];
zonal_mask = [zonal_mask zeros(5,3)];
zonal_mask = [zonal_mask;zeros(3,8)];

threshold_mask = [1 1 0 1 1 0 0 0;1 1 1 1 0 0 0 0;1 1 0 0 0 0 0 0;1 0 0 0 0 0 0 0];
threshold_mask = [threshold_mask;zeros(1,8);0 1 0 0 0 0 0 0;zeros(2,8)];

figure();
imshow(haar(image,2),[]);

function image_p = get_part(image,part)
    [row,col] = size(image);
    if part == 1
        image_p = image(1:row/2,1:col/2);
    elseif part == 2
        image_p = image(1:row/2,col/2+1:col);
    elseif part == 3
        image_p = image(row/2+1:row,1:col/2);
    elseif part == 4
        image_p = image(row/2+1:row,col/2+1:col);
    end
    
end
function image_ = haar(image,round)
    [row,col] = size(image);
    image_ = image;
    
    for r = 1:round
        image_(1:row/(2^r),1:col/(2^r)) = (get_part(image,1)+get_part(image,2)+get_part(image,3)+get_part(image,4))/2;
        image_(1:row/(2^r),col/(2^r)+1:col/(2^(r-1))) = (get_part(image,1)-get_part(image,2)+get_part(image,3)-get_part(image,4))/2;
        image_(row/(2^r)+1:row/(2^(r-1)),1:col/(2^r)) = (get_part(image,1)+get_part(image,2)-get_part(image,3)-get_part(image,4))/2;
        image_(row/(2^r)+1:row/(2^(r-1)),col/(2^r)+1:col/(2^(r-1))) = (get_part(image,1)-get_part(image,2)-get_part(image,3)+get_part(image,4))/2;
        image = image_(1:row/(2^r),1:col/(2^r));
    end
    image_ = image_(1:row/(2^round),1:col/(2^round));
end


function m = get_matrix(n)
    m = zeros(n);
    for u = 0:n-1
        for x = 0:n-1
            m(u+1,x+1)=alpha(u,n)*cos(((2*x+1)*u*pi)/(2*n));
        end
    end
end

function image_ = DCT(image,size_,mask)
    image = convert(image,size_);
    [row,col]=size(image);
    image_ = im2double(zeros(row,col));
    
    m = get_matrix(size_);
    
    for r = 1:size_:row+1-size_
        for c = 1:size_:col+1-size_
            image_(r:r+size_-1,c:c+size_-1) = m*image(r:r+size_-1,c:c+size_-1)*m';
        end
    end
    
    for r = 1:size_:row+1-size_
        for c=1:size_:col+1-size_
            image_(r:r+size_-1,c:c+size_-1) = image_(r:r+size_-1,c:c+size_-1).*mask;
        end
    end
end
function image_ = IDCT(image,size_)
    [row,col]=size(image);
    image_ = im2double(zeros(row,col));
    m = get_matrix(size_);
    for r = 1:size_:row+1-size_
        for c = 1:size_:col+1-size_
            image_(r:r+size_-1,c:c+size_-1) = m'*image(r:r+size_-1,c:c+size_-1)*m;
        end
    end
end
function v = alpha(u,size_)
    if u == 0
        v = sqrt(1/size_);
    else
        v = sqrt(2/size_);
    end
end
function image_ = convert(image,size_)
    [row,col]=size(image);
    size_x = mod(row,size_);
    size_y = mod(col,size_);
    image_ = [image zeros(row,size_y)];
    image_ = [image_;zeros(size_x,col+size_y)];
    image_ = im2double(image_);
end