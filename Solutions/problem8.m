%{
    Problem 8
    by Fanyong Xue
	Student ID:515030910443
    Morpholigical Processing
    Note: this problem may take too some time to finish
%}

%% Main Part
clear all
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/Fig0929(a)(text_image).tif');
figure();
imshow(image);

b=erode(image,51,1);
figure();
imshow(b);
c = dilate(b,51,1);
figure();
imshow(c);

d = rd(b,image,8);
figure();
imshow(d);

e=filling_holes(image);
figure();
imshow(e);

f = border_clearing(image);
figure();
imshow(f);

%% Function Part
function image_ = dilate(image,mask_size_x,mask_size_y)
    [row,col] = size(image);
    %{
    image_ = zeros(row,col);
    
    image = [zeros(row,floor(mask_size_y/2)) image zeros(row,floor(mask_size_y/2))];
    image = [zeros(floor(mask_size_x/2),col+mask_size_y-1);image;zeros(floor(mask_size_x/2),col+mask_size_y-1)];

    for r = 1:row
        for c = 1:col
            image_(r,c) = logical(max(max(image(r:r+mask_size_x-1,c:c+mask_size_y-1))));
        end
    end
    %}
    image_ = zeros(row+mask_size_x-1,col+mask_size_y-1);%append
    for r = 1:row
        for c = 1:col
            if image(r,c)
                image_(r:r+mask_size_x-1,c:c+mask_size_y-1) = ones(mask_size_x,mask_size_y);
            end
        end
    end
    image_ = logical(image_(ceil(mask_size_x/2):row+floor(mask_size_x/2),ceil(mask_size_y/2):col+floor(mask_size_y/2)));
end

function image_ = erode(image,mask_size_x,mask_size_y)
    [row,col] = size(image);
    image_ = zeros(row,col);
    
    image = [ones(row,floor(mask_size_y/2)) image ones(row,floor(mask_size_y/2))];
    image = [ones(floor(mask_size_x/2),col+mask_size_y-1);image;ones(floor(mask_size_x/2),col+mask_size_y-1)];

    for r = 1:row
        for c = 1:col
            image_(r,c) = logical(min(min(image(r:r+mask_size_x-1,c:c+mask_size_y-1))));
        end
    end
end

function image_ = rd(F,G,f)
    figure(f);
    pre_DGF = logical(and_image(dilate(F,3,3),G));
    DGF = logical(and_image(dilate(pre_DGF,3,3),G));
    while image_different(pre_DGF,DGF)
        pre_DGF = DGF;
        DGF = logical(and_image(dilate(pre_DGF,3,3),G));
        imshow(DGF);
    end
    image_ = DGF;
    close(figure(f));
end

function image_ = re(F,G)
    pre_DGF = logical(or_image(erode(F,3,3),G));
    DGF = logical(or_image(erode(pre_DGF,3,3),G));
    while image_different(pre_DGF,DGF)
        pre_DGF = DGF;
        DGF = logical(and_image(erode(pre_DGF,3,3),G));
    end
    image_ = DGF;
end

function image_ = and_image(image1,image2)
    image_ = image1&image2;
end
function image_ = or_image(image1,image2)
    image_ = image1|image2;
end
function v = image_different(image1,image2)
    v = sum(sum(image1-image2));
end
function image_ = filling_holes(image)
    image_c = ~image;
    [row,col] = size(image);
    image_ = ones(row,col) - image;
    image_(2:row-1,2:col-1) = zeros(row-2,col-2);
    image_ = logical(image_);
    image_h = rd(image_,image_c,9);
    image_ = ~image_h;
    
end
function image_ = border_clearing(image)
    [row,col] = size(image);
    image_f = image;
    image_f(2:row-1,2:col-1) = zeros(row-2,col-2);
    image_f = logical(image_f);
    
    image_ = image - rd(image_f,image,7);
    
end