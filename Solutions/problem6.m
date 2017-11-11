%{
    Problem 6
    by Fanyong Xue 
    Student ID:515030910443
    Geometric transform
%}

%% Main Part
clear all
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/ray_trace_bottle.tif');
figure('name','rotate');
subplot(311);
imshow(image);
title('(a) Original Image');
subplot(312);
imshow(rotate(image,pi/6,'nearest neighbor'));
title('Rotate (a) by pi/6 using nearest neighbor');
subplot(313);
imshow(rotate(image,pi/6,'bilinear'));
title('Rotate (a) by pi/6 using bilinear');


figure('name','translate');
subplot(311);
imshow(image);
title('(a) Original Image');
subplot(312);
imshow(translate(image,500,500));
title('Translate (a) by x:500 and y:500');
subplot(313);
imshow(translate(image,-500,-500));
title('Translate (a) by x:-500 and y:-500');


figure('name','scale');
subplot(221);
imshow(image);
title('(a) Original Image');
subplot(222);
image_b = scale(image,0.25,0.25,'bilinear');
imshow(image_b);
title('(b) Scale (a) to 0.25*0.25');
subplot(223);
imshow(scale(image_b,2,2,'nearest neighbor'));
title('Scale (b) to 2*2 using nearest neighbor');
subplot(224);
imshow(scale(image_b,2,2,'bilinear'));
title('Scale (b) to 2*2 using bilinear');

%% Function Part

% rotate the point(v,w) angle degrees to point(x,y)
function [x,y]=rotate_axis(v,w,angle)
    x = round(v*cos(angle)-w*sin(angle));
    y = round(v*sin(angle)+w*cos(angle));
end
%{ 
    totate the image angle degrees with specified interpolation
    interpolation:'nearest neighbor' or 'bilinear'
%}
function image_ = rotate(image,angle,interpolation)
    [row,col] = size(image);
    x1 = 0;
    y1 = 0;
    [x2,y2] = rotate_axis(0,col-1,angle);
    [x3,y3] = rotate_axis(row-1,col-1,angle);
    [x4,y4] = rotate_axis(row-1,0,angle);
    y_r = max([y1 y2 y3 y4]) - min([y1 y2 y3 y4]);
    x_r = max([x1 x2 x3 x4]) - min([x1 x2 x3 x4]);
    image_ = -1*(ones(x_r+1,y_r+1));
    
    x_min = min([x1 x2 x3 x4]);
    y_min = min([y1 y2 y3 y4]);
    for x = 1:row
        for y = 1:col
            [x_,y_] = rotate_axis(x-1,y-1,angle);
            image_(x_-x_min+1,y_-y_min+1) = image(x,y);
        end
    end
    
    edge = zeros(x_r+1,2);
    for x = 1:x_r+1
        for y = 1:y_r+1
            if image_(x,y) == -1
                image_(x,y) = 0;
                continue;
            else
                break;
            end
        end
        edge(x,1) = y;
        for y = y_r+1:-1:1
            if image_(x,y) == -1
                image_(x,y) = 0;
            else
                break;
            end
        end
        edge(x,2) = y;
    end
    if strcmp(interpolation,'nearest neighbor')
        for x = 1:x_r+1
            for y = 1:y_r+1
                if image_(x,y) == -1 
                    [x_ori,y_ori] = rotate_axis(x-1+x_min,y-1+y_min,-1*angle);
                    v = find_nearest(image,x_ori+1,y_ori+1);
                    image_(x,y) = v(1);
                end
            end
        end
    elseif strcmp(interpolation,'bilinear')
        for x = 1:x_r+1
            for y = 1:y_r+1
                if image_(x,y) == -1 
                    [x_ori,y_ori] = rotate_axis(x-1+x_min,y-1+y_min,-1*angle);
                    v = find_nearest(image,x_ori+1,y_ori+1);
                    image_(x,y) = mean(v);
                end
            end
        end
    end
    
    image_ = uint8(image_);
end
% find the nearest point(up,down,left,right if have) which value is not -1
% and store them on v
function v = find_nearest(image,x,y)
    [row,col] = size(image);
    v = 0;
    found = 0;
    %up
    if x > 1 && image(x-1,y)~=-1
        v = image(x-1,y);
        found = 1;
    end
    %down
    if x<row && image(x+1,y)~=-1
        if found
            v = [v image(x+1,y)];
        else
            v = image(x+1,y);
            found =1;
        end
        
    end
    %left
    if y<col && image(x,y+1)~=-1
        if found
            v = [v image(x,y+1)];
        else
            v = image(x,y+1);
            found =1;
        end
    end
    %right
    if y>1 && image(x,y-1)~=-1
        v = [v image(x,y-1)];
        if found
            v = [v image(x,y-1)];
        else
            v = image(x,y-1);
        end
    end
end

% translate the image by tx in x dimension and ty in y dimension

function image_ = translate(image,tx,ty)
    [row,col] = size(image);
    
    if tx>0
        image_ = [zeros(row,tx) image];
    else
        image_ = [image zeros(row,-1*tx)];
    end
    
    if ty>0
        image_ = [zeros(ty,col+abs(tx));image_];
    else
        image_ = [image_;zeros(abs(ty),col+abs(tx))];
    end
end
%{ 
    scale the image by cx in x dimension and cy in y dimension
    with specified interpolation
    interpolation:'nearest neighbor' or 'bilinear'
%}
function image_ = scale(image,cx,cy,interpolation)
    [row,col] = size(image);
    row_ = ceil(row*cx);
    col_ = ceil(col*cy);
    image_ = -1*ones(row_,col_);
    
    for x = 1:row
        for y = 1:col
            image_(ceil(x*cx),ceil(y*cy)) = image(x,y);
        end
    end
    if cx>1 || cy>1
        if strcmp(interpolation,'nearest neighbor')
            for x = 1:row_
                for y = 1:col_
                    if image_(x,y)==-1
                        x_ori = ceil(x/cx);
                        y_ori = ceil(y/cy);
                        v = find_nearest(image,x_ori,y_ori);
                        image_(x,y) = v(1);
                    end
                end
            end
        elseif strcmp(interpolation,'bilinear')
            for x = 1:row_
                for y = 1:col_
                    if image_(x,y)==-1
                        x_ori = ceil(x/cx);
                        y_ori = ceil(y/cy);
                        v = find_nearest(image,x_ori,y_ori);
                        image_(x,y) = mean(v);
                    end
                end
            end
        end
    end
    image_ = uint8(image_);
end