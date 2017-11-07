%{
%}

image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/ray_trace_bottle.tif');
image = im2double(image);




function [x,y]=rotate_axis(v,w,angle)
    x = round(v*cos(angle)-w*sin(angle));
    y = round(v*sin(angle)+w*cos(angle));
end

function image_ = rotate(image,anlge)
    [row,col] = size(image);
    x1 = 0;
    y1 = 0;
    [x2,y2] = rotate_axis(0,col-1,angle);
    [x3,y3] = rotate_axis(row-1,col-1,angle);
    [x4,y4] = rotate_axis(row-1,0,angle);
    y_r = max([y1 y2 y3 y4]) - min([y1 y2 y3 y4]);
    x_r = max([x1 x2 x3 x4]) - min([x1 x2 x3 x4]);
    image_ = uint8(zeros(x_r,y_r));
    
    x1_ = 0;
    y1_ = 0;
    for x = 1:row
        for y = 1:col
            [x_,y_] = rotate_axis(x-1,y-1,angle);
            image_(x_+x1_,y_+y1_) = image(x,y);
        end
    end
end