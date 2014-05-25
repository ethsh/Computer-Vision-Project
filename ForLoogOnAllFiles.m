close all;
clear parameters;

for i = 1:25
    if (i < 10)
        num = ['0' num2str(i)];
    else
        num = num2str(i);
    end
    
    [im, map] = imread(['balls/MVC-0' num 'F.JPG']);
    im_hsv = rgb2hsv(im);
    
%     figure; imshow(im_hsv(:,:,2));
    figure; imshow(im);
end