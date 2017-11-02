image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/skeleton_orig.tif');
[row,col] = size(image);
c_image = image;
d_image = image;
e_image = image;
f_image = image;
g_image = image;
h_image = image;
double_image = double(image);
lap_image = zeros(row,col);
mask = [-1 -1 -1;-1 8 -1;-1 -1 -1];

%append_image = [zeros(row,1) double_image zeros(row,1)];
%append_image = [zeros(1,col+2); append_image ;zeros(1,col+2)];

%d
x_mask = [-1 -2 -1;0 0 0;1 2 1];
y_mask = [-1 0 1;-2 0 2;-1 0 1];


for r = 2:row-1
    for c = 2:col-1
        lap_image(r,c) = image(r,c)+sum(sum(double_image(r-1:r+1,c-1:c+1).*mask));
        c_image(r,c) = image(r,c) - lap_image(r,c);
        d_image(r,c) = abs(sum(sum(double_image(r-1:r+1,c-1:c+1).*x_mask)))+abs(sum(sum(double_image(r-1:r+1,c-1:c+1).*y_mask)));
    end
end



min = lap_image(1,1);
max = lap_image(1,1);

for r = 1:row
    for c = 1:col
        if min > lap_image(r,c)
            min = lap_image(r,c);
        end
        if max < lap_image(r,c)
            max = lap_image(r,c);
        end
    end
end

max = max - min;
lap_image = lap_image - min;
lap_image = (lap_image/max);
lap_image = lap_image+0.8;

%d_image

for r = 3:row-2
    for c = 3:col-2
        e_image(r,c) = mean(mean(double_image(r-2:r+2,c-2:c+2)));
    end
end

%f_image g_image h_image


f_image = e_image.*c_image;
g_image = f_image+image;
h_image = sqrt(double(g_image));


figure(1);
subplot(241);
imshow(double(zeros(row,col))+0.8);
subplot(242);
imshow(lap_image);

subplot(245);
imshow(c_image);
subplot(246);
imshow(d_image);

subplot(243);
imshow(e_image);

subplot(244);
imshow(f_image);

subplot(247);
imshow(g_image);

subplot(248);
imshow(h_image);





%{
function image_f = fourier_transform(image,rows,cols)
    [row,col]=size(image);
    
    image = [image zeros(row,cols-col)];
    image = [image;zeros(rows-row,cols)];
    
    image_f = double(zeros(rows,cols));
    base_x_matrix = zeros(rows,1);
    base_y_matrix = zeros(cols,1);
    
    
    for x = 1:rows
        base_x_matrix(x) = exp(1i*2*pi*(x-1)/rows);
    end
    
    for y = 1:cols
        base_y_matrix(y) = exp(1i*2*pi*(y-1)/cols);
    end
    
    
    image_f(1,1) = sum(sum(image));
    
    x_sum = sum(double(image'));
    y_sum = sum(double(image));
    
    for r = 2:rows
        x_sum = x_sum .* base_x_matrix';
        image_f(r,1) = sum(x_sum);
    end
        
    for c = 2:cols
        y_sum = y_sum .* base_y_matrix';
        image_f(1,c) = sum(y_sum);
    end
    image_f(1,56)
    image_f(56,1)
    end
%}
