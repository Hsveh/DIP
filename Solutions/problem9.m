%{
    Problem 9
    by Fanyong Xue
	Student ID:515030910443
    Image segmentation
%}
%% Main Part
clear all
image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/building.tif');
image2 = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/polymersomes.tif');

Sobel_x = double([-1 -2 -1;0 0 0;1 2 1]);
Sobel_y = double([-1 0 1;-2 0 2;-1 0 1]);
Prewitt_x = double([-1 -1 -1;0 0 0;1 1 1]);
Prewitt_y = double([-1 0 1;-1 0 1;-1 0 1]);
Roberts_x = double([-1 0;0 1]);
Roberts_y = double([0 -1;1 0]);


[k1,k2]=detectors(image,Sobel_x,Sobel_y);
figure('name','Sobel x');
imshow(k1,[]);
figure('name','Sobel y');
imshow(k2,[]);
figure('name','Sobel');
imshow(k1+k2,[]);


[k1_,k2_]=detectors(image,Prewitt_x,Prewitt_y);
figure('name','Prewitt x');
imshow(k1_,[]);
figure('name','Prewitt y');
imshow(k2_,[]);
figure('name','Prewitt');
imshow(k1_+k2_,[]);


[k1__,k2__]=detectors(image,Roberts_x,Roberts_y);
figure('name','Roberts x');
imshow(k1__,[]);
figure('name','Roberts y');
imshow(k2__,[]);
figure('name','Roberts');
imshow(k1__+k2__,[]);

figure('name','Marr Hildreth');
image_m1=Marr_Hildreth(image,4,25);
imshow(image_m1);

figure('name','Canny');
imshow(Canny(image));
figure('name','Basic Global Thresholding');
imshow(Basic_Global_Thresholding(image2,255));
figure('name','Otsu');
imshow(Otsu(image2));


%% Function Part

function [image_x,image_y] = detectors(image,mask_x,mask_y)
    
    image = double(image);
    image_x = image;
    image_y = image;
    [row,col] = size(image);
    [m_x,m_y] = size(mask_x);
    if m_x <= 2 && m_y <=2
        mask_x = [zeros(m_x,1) mask_x];
        mask_x = [zeros(1,1+m_y);mask_x];
        mask_y = [zeros(m_x,1) mask_y];
        mask_y = [zeros(1,1+m_y);mask_y];
    end
    [m_x,m_y] = size(mask_x);
    append_x = floor(m_x/2);
    append_y = floor(m_y/2);
    %image append
    image = [ones(append_y,1)*image(1,:);image;ones(append_y,1)*image(end,:)];
    image = [image(:,1)*ones(1,append_x) image image(:,end)*ones(1,append_x)];
    
    for r = 1:row
        for c = 1:col
            image_x(r,c) = abs(sum(sum(image(r:r+m_x-1,c:c+m_y-1).*mask_x))); 
        end
    end

    for r = 1:row
        for c = 1:col
            image_y(r,c) = abs(sum(sum(image(r:r+m_x-1,c:c+m_y-1).*mask_y))); 
        end
    end

end

function image__ = Marr_Hildreth(image,sigma,n)
    g = zeros(n);
    image = double(image);
    image_=image;
    [row,col] = size(image);
    image__= zeros(row,col);
    for r = -1*floor(n/2):floor(n/2)
        for c = -1*floor(n/2):floor(n/2)
            g(r+floor(n/2)+1,c+floor(n/2)+1) = ((r*r+c*c-2*sigma*sigma)/(sigma^4))*exp(-1*(r*r+c*c)/(2*sigma*sigma));
        end
    end
    
    image = [ones(floor(n/2),1)*image(1,:);image;ones(floor(n/2),1)*image(end,:)];
    image = [image(:,1)*ones(1,floor(n/2)) image image(:,end)*ones(1,floor(n/2))];
    
    for r = 1:row
        for c = 1:col
            image_(r,c) = sum(sum(image(r:r+n-1,c:c+n-1).*g));
        end
    end
    t = 0.04*max(max(image_));
    %zero
    
    for r = 2:row-1
        
        for c = 2:col-1
            zeros_ = 0;
            if po_ne(image_(r+1,c))+po_ne(image_(r-1,c))==0 && abs(image_(r+1,c)-image_(r-1,c))>=t
                zeros_ = zeros_+1;
            end
            if po_ne(image_(r,c+1))+po_ne(image_(r,c-1))==0 && abs(image_(r,c+1)-image_(r,c-1))>=t
                zeros_ = zeros_+1;
            end
            if po_ne(image_(r+1,c+1))+po_ne(image_(r-1,c-1))==0 && abs(image_(r+1,c+1)-image_(r-1,c-1))>=t
                zeros_ = zeros_+1;
            end
            if po_ne(image_(r+1,c-1))+po_ne(image_(r-1,c+1))==0 && abs(image_(r+1,c-1)-image_(r-1,c+1))>=t
                zeros_ = zeros_+1;
            end
            if zeros_>=1
                image__(r,c)=255;
            else
                image__(r,c) = 0;
            end
            
        end
    end
end

function v_ = po_ne(v)
    if v > 0
        v_ = 1;
    elseif v<0
        v_ = -1;
    else
        v_ =2;
    end
end

function image_res = Canny(image)
    n = 25;
    sigma=4;
    g = zeros(n);
    image = double(image);
    image_=image;
    [row,col] = size(image);
    for r = -1*floor(n/2):floor(n/2)
        for c = -1*floor(n/2):floor(n/2)
            g(r+floor(n/2)+1,c+floor(n/2)+1) = exp(-1*(r*r+c*c)/(2*sigma*sigma));
        end
    end
    image = [ones(floor(n/2),1)*image(1,:);image;ones(floor(n/2),1)*image(end,:)];
    image = [image(:,1)*ones(1,floor(n/2)) image image(:,end)*ones(1,floor(n/2))];
    
    for r = 1:row
        for c = 1:col
            image_(r,c) = sum(sum(image(r:r+n-1,c:c+n-1).*g));
        end
    end
    Roberts_x = double([-1 0;0 1]);
    Roberts_y = double([0 -1;1 0]);
    [gx,gy] = detectors(image_,Roberts_x,Roberts_y);
    M = sqrt(gx.^2+gy.^2);
    alpha = atan(gy./gx);
    gn = M;
    for r = 2:row-1
        for c = 2:col-1
            p = part(alpha(r,c));
            if p == 1
                nei1 = M(r-1,c);
                nei2 = M(r+1,c);
            elseif p==2
                nei1 = M(r-1,c-1);
                nei2 = M(r+1,c+1);
            elseif p==3
                nei1 = M(r,c+1);
                nei2 = M(r,c-1);
            else
                nei1 = M(r-1,c+1);
                nei2 = M(r+1,c-1);
            end
            if M(r,c)<nei1 || M(r,c)<nei2
                gn(r,c)=0;
            else
                gn(r,c)=M(r,c);
            end
        end
    end
    %gn=uint8(gn);
    tl=185;
    th=370;
    gh=gn*0;
    gl=gn*0;
    for r = 1:row
        for c = 1:col
            if gn(r,c)>=th
                gh(r,c)=gn(r,c);
                gl(r,c)=gn(r,c);
            elseif gn(r,c)>=tl
                gl(r,c)=gn(r,c);
            end
        end
    end
    image_res = zeros(row,col);
    go = logical(zeros(row,col));
    for r = 1:row
        for c = 1:col
            if gh(r,c)>0 && image_res(r,c)==0
                image_res(r,c) = gh(r,c);
                [image_res,go]=go_through(image_res,go,gl,r,c,row,col);
                
            end
        end
    end
    
end

function [image_,go] = go_through(image_,go,image,x,y,row,col)
    if x>1&&image(x-1,y)>0&&~go(x-1,y)
        image_(x-1,y)=image(x-1,y);
        go(x-1,y)=1;
        [image_,go]=go_through(image_,go,image,x-1,y,row,col);
    end
    if y>1&&image(x,y-1)>0&&~go(x,y-1)
        image_(x,y-1)=image(x,y-1);
        go(x,y-1)=1;
        [image_,go]=go_through(image_,go,image,x,y-1,row,col);
       
    end
    
    if x<row&&image(x+1,y)>0&&~go(x+1,y)
        image_(x+1,y)=image(x+1,y);
        go(x+1,y)=1;
        [image_,go]=go_through(image_,go,image,x+1,y,row,col);
       
    end
    if y<col&&image(x,y+1)>0&&~go(x,y+1)
        image_(x,y+1)=image(x,y+1);
        go(x,y+1)=1;
        [image_,go]=go_through(image_,go,image,x,y+1,row,col);
       
    end
    if x>1&&y>1&&image(x-1,y-1)>0&&~go(x-1,y-1)
        image_(x-1,y-1)=image(x-1,y-1);
        go(x-1,y-1)=1;
        [image_,go]=go_through(image_,go,image,x-1,y-1,row,col);
       
    end
    if x<row&&y<col&&image(x+1,y+1)>0&&~go(x+1,y+1)
        image_(x+1,y+1)=image(x+1,y+1);
        go(x+1,y+1)=1;
        [image_,go]=go_through(image_,go,image,x+1,y+1,row,col);
       
    end
    if x>1&&y<col&&image(x-1,y+1)>0&&~go(x-1,y+1)
        image_(x-1,y+1)=image(x-1,y+1);
        go(x-1,y+1)=1;
        [image_,go]=go_through(image_,go,image,x-1,y+1,row,col);
       
    end
    if x<row&&y>1&&image(x+1,y-1)>0&&~go(x+1,y-1)
        image_(x+1,y-1)=image(x+1,y-1);
        go(x+1,y-1)=1;
        [image_,go]=go_through(image_,go,image,x+1,y-1,row,col);
       
    end
end
function v = part(angle)
    if (-pi/2 <= angle && angle <= -3*pi/8)||(pi/2 >= angle && angle >= 3*pi/8)
        v= 1;
    elseif -3*pi/8<= angle && angle <=-1*pi/8
        v=2;
    elseif -1*pi/8<=angle && angle<=1*pi/8
        v = 3;
    else
        v=4;
    end
end

function image_ = Basic_Global_Thresholding(image,delta)
    [row,col] = size(image);
    g1=zeros(row,col);
    g2=zeros(row,col);
    m = mean(mean(image));
    prevm = 0;
    while abs(m-prevm)>delta
        for r = 1:row
            for c = 1:col
                if image(r,c)>m
                    g1(r,c) = image(r,c);
                else
                    g2(r,c) = image(r,c);
                end
            end
        end
        prevm = m;
        m = (mean(mean(g1))+mean(mean(g2)))/2;
    end
    image_ = zeros(row,col);
    for r = 1:row
        for c = 1:col
            if image(r,c)>m
                image_(r,c)=image(r,c);
            end
        end
    end
end
function image_=Otsu(image)
    
    [row,col] = size(image);
    
    image_ = zeros(row,col);
    p = get_histogram(image);
    p = double(p);
    p = p./(row*col);
    P1 = zeros(256,1);
    P1(1) = p(1);
    for k=2:256
        P1(k) = P1(k-1)+p(k);
    end

    m = zeros(256,1);
    m(1)=0;
   
    for k=2:256
        m(k) = m(k-1)+(k-1)*p(k);
    end
    mg = m(256);
    max_v = 0;
    max_k=zeros(256,1);
    max_k(1) = 1;
    num_k=1;
    for k = 1:256
        sigma = (mg*P1(k)-m(k))^2/(P1(k)*(1-P1(k)));
        if sigma > max_v
            max_k(1)=k;
            num_k=1;
            max_v = sigma;
        elseif sigma == max_v
            num_k = num_k+1;
            max_k(num_k) = k;
        end
    end
    max_k_ = mean(max_k(1:num_k));
    for r = 1:row
        for c = 1:col
            if image(r,c)>max_k_
                image_(r,c)=image(r,c);
            end
        end
    end
end

function histogram = get_histogram(image)
    histogram = zeros(256,1);
    [row,col]=size(image);
    for r = 1:row
        for c = 1:col
            gray = image(r,c);
            histogram(gray+1)=histogram(gray+1)+1;
        end
    end
end