% Problem 1
% by Fanyong Xue
% Student ID:515030910443
% Histogram Equalizatio

%% Main Part
clear all
image1 = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/Fig1.jpg');
image2 = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/Fig2.jpg');

[histogram1,histogram_e1,transfer_f1,image_e1] = histogram_equalization(image1);
[histogram2,histogram_e2,transfer_f2,image_e2] = histogram_equalization(image2);

plot_data(image1,image_e1,histogram1,histogram_e1,transfer_f1);
plot_data(image2,image_e2,histogram2,histogram_e2,transfer_f2);

%% Functions Part

% get histogram of image
% image: get histogram of it
% histogram: the histogram of image
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

%do the histogram_equalization for image
%image: do histogram_equalization for it
%histogram: original histogram; histogram_e: histogram after histogram
%equalizatio; transfer_f: transfer function; image_e: image after histogram
%equalizatio
function [histogram,histogram_e,transfer_f,image_e] = histogram_equalization(image)
    [row,col]=size(image);
    transfer_f = zeros(256,1);
    histogram = get_histogram(image);
    transfer_f(1) = 256*histogram(1)/(row*col);
    
    for i = 2:256
        transfer_f(i) = transfer_f(i-1)+255*histogram(i)/(row*col);
    end
    transfer_f = round(transfer_f);
    
    image_e = image;
    for r = 1:row
        for c = 1:col
            image_e(r,c)=transfer_f(image(r,c)+1);
        end
    end
    histogram_e = get_histogram(image_e);
end

%plot data
%image:original image; image_e: image after histogram equalizatio; histogram: original histogram; 
%histogram_e: histogram after histogram equalizatio; transfer_f: transfer function
function plot_data(image,image_e,histogram,histogram_e,transfer_f)
    figure();
    subplot(2,3,1);
    imshow(image);
    title("Original Image");
    subplot(2,3,2);
    imshow(image_e);
    title("Image(Histogram Equalization)");
    subplot(2,3,3);
    bar(histogram);
    title("Histogram");
    subplot(2,3,4);
    bar(histogram_e);
    title("Histogram(Equalization)");
    subplot(2,3,5);
    plot(transfer_f);
    title("Transfer Funciton");
end
