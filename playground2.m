clc; close all; clear all;

%% Histogram Calculation - through Ball001.JPG !
[im, map] = imread('balls/MVC-001F.JPG');

im_hsv = rgb2hsv(im);
im_hs = im_hsv(:,:,1:2);

im_hs = im_hs*255;

ball_center = [350 240];
ball_winsize = 40*[1 1];

Ball = impart(im_hs,ball_center,ball_winsize);

bins = 32;
binsize = 256 / bins;

bin1 = ceil(double(Ball(:,:,1))/(binsize+1))+1;
bin2 = ceil(double(Ball(:,:,2))/(binsize+1))+1;
% bin3 = ceil(double(Ball(:,:,3))/(binsize+1))+1;

% H = zeros(bins,bins,bins);
H = zeros(bins);
for i=1:size(Ball,1)
    for j=1:size(Ball,2)
%         H(bin1(i,j),bin2(i,j),bin3(i,j)) = H(bin1(i,j),bin2(i,j),bin3(i,j)) + 1;
          H(bin1(i,j),bin2(i,j)) = H(bin1(i,j),bin2(i,j)) + 1;
%           H(bin1(i,j)) = H(bin1(i,j)) + 1;
    end
end

%% Histogram training through MVC009 too.

[im2, map2] = imread('balls/MVC-009F.JPG');

im_hsv2 = rgb2hsv(im2);
im_hs2 = im_hsv2(:,:,1:2);

im_hs2 = im_hs2*255;

ball_center2 = [434 285];
ball_winsize2 = 30 * [1 1];

Ball2 = impart(im_hs2,ball_center2,ball_winsize2);
bin1 = ceil(double(Ball2(:,:,1))/(binsize+1))+1;
bin2 = ceil(double(Ball2(:,:,2))/(binsize+1))+1;
% bin3 = ceil(double(Ball(:,:,3))/(binsize+1))+1;

for i=1:size(Ball2,1)
    for j=1:size(Ball2,2)
%         H(bin1(i,j),bin2(i,j),bin3(i,j)) = H(bin1(i,j),bin2(i,j),bin3(i,j)) + 1;
          H(bin1(i,j),bin2(i,j)) = H(bin1(i,j),bin2(i,j)) + 1;
%           H(bin1(i,j)) = H(bin1(i,j)) + 1;
    end
end

H = H ./ max(H(:)) * 255;

%% Read new img
[img,map_img] = imread('balls/MVC-025F.JPG');

figure; imshow(img);

img_hsv = rgb2hsv(img);
img_hs = img_hsv(:,:,1:2);

img_hs = img_hs*255;

%% Backprojecting on a new img:
new_img = BackProjection(img_hs,H,bins,2);
figure; imshow(new_img,map_img);

%% Thresholding the image

thresh = 10;

thresh_img = zeros(size(new_img));
thresh_img(new_img>thresh)=new_img(new_img>thresh);
figure; imshow(thresh_img, map_img)

%% Full circle search
% r = 10:10:100;
r = 40;
for k = 1:length(r)
    r_t = r(k);
    filt = zeros(2*r_t);
    theta = 0:0.05:2*pi;
    r2 = 0:0.05:r_t;
    for i=1:length(theta)
        for j = 1:length(r2)
            filt(ceil(r_t+r2(j)*sin(theta(i))),ceil(r_t+r2(j)*cos(theta(i)))) = 1/(pi*r_t^2);
        end
    end
    
    res = filter2(filt, thresh_img);
    res2 = filter2(filt, new_img);
    
    
    figure; imshow(res2, map);
%     figure; imshow(peaks2, map);
    [m, ~] = max(res2(:));
    title(['r = ' num2str(r_t) ', NOT Thresholded image, the max is ' num2str(m)]);
    
    figure; imshow(res, map);
%     figure; imshow(peaks, map);
    [m, ~] = max(res(:));
    title(['r = ' num2str(r_t) ', Thresholded image, the max is ' num2str(m)]);

    
    
    
%     [peaks, locs] = findpeaks(res, 'MINPEAKDISTANCE', r_t);
    BW = imregionalmax(res);
    [m, i] = max(res(:));
    [y, x] = ind2sub(size(res), i);
    
    %% find peaksun-efficiently
%     peaks = zeros(size(res));
%     peaks2 = zeros(size(res2));
%     
%     for i = 2:size(res, 1)-1
%         for l =2:size(res, 2)-1
%             if (res(i,l)>res(i-1,l) && res(i,l)>res(i+1,l) && res(i,l)>res(i,l-1) && res(i,l)>res(i,l+1) && ...
%                 res(i,l)>res(i-1,l-1) && res(i,l)>res(i-1,l+1) && res(i,l)>res(i+1,l-1) && res(i,l)>res(i+1,l+1))
%                 peaks(i,l) = 1;
%             end
%             
%             if (res2(i,l)>res2(i-1,l) && res2(i,l)>res2(i+1,l) && res2(i,l)>res2(i,l-1) && res2(i,l)>res2(i,l+1) && ...
%                 res2(i,l)>res2(i-1,l-1) && res2(i,l)>res2(i-1,l+1) && res2(i,l)>res2(i+1,l-1) && res2(i,l)>res2(i+1,l+1))
%                 peaks2(i,l) = 1;
%             end
%         end
%     end
%     
    
%     
%     peaks_vals = res(peaks==1);
%     peaks_vals2 = res2(peaks2==1);
%     
%     [m, i] = max(peaks_vals(:));
%     [y, x] = ind2sub(size(peaks_vals), i)
%     peaks_vals(y, x) = 0;
%     [m, i] = max(peaks_vals(:));
%     [y, x] = ind2sub(size(peaks_vals), i)
%     
%     [m2, i2] = max(peaks_vals2(:));
%     [y2, x2] = ind2sub(size(peaks_vals2), i2)
%     peaks_vals2(y2, x2) = 0;
%     [m, i] = max(peaks_vals(:));
%     [y2, x2] = ind2sub(size(peaks_vals), i2)
    
end

%% Filter
sigma = 10;
h = fspecial('gaussian', 15, sigma);
filtered_img = imfilter(new_img, h);
figure; imshow(filtered_img, map_img);

%% Circle finding

edged_img = edge(filtered_img,'canny');
figure; imshow(edged_img, map_img);

radiusRange = [20 40];
%[centers,radii] = imfindcircles(filtered_img,radiusRange);
[centers,radii] = imfindcircles(edged_img,radiusRange, 'Sensitivity',0.9);

%% Circle plotting
figure; imshow(img);
hold on;
viscircles(centers, radii);