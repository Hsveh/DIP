image = imread('/Users/xuefanyong/Documents/GitHub/DIP/Solutions/images/Fig0503.tif');
[row,col] = size(image);

function n = uniform_noise(row,col,low,high)
    %n = zeros(row,col);
    
    n = uint8(low + (high-low)*rand([row col]));
    
end

function n = gaussian_noise(row,col,ex,va)
    n_ = uniform_noise(row,col,0,255);
    n = uint8(zeros(row,col));
    n_ = double(n_);
    n_ = n_./255;
    
    for r = 1:row
        for c = 1:col
            n(r,c) = sqrt(-2*log(n_(r,c))*cos(2*pi*n_(r,mod(c+1,col)+1)));
        end
    end
    
    n = n*va;
    n = n+ex;
end

