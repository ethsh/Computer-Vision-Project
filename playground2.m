clc; close all; clear all;

%% Histogram Calculation - through Ball001.JPG !
[im, map] = imread('ball01.JPG');

ball_center = [350 240];
ball_winsize = 40*[1 1];

Ball = impart(im,ball_center,ball_winsize);

bins = 20;
binsize = 256 / bins;

bin1 = ceil(double(Ball(:,:,1))/(binsize+1))+1;
bin2 = ceil(double(Ball(:,:,2))/(binsize+1))+1;
bin3 = ceil(double(Ball(:,:,3))/(binsize+1))+1;

H = zeros(bins,bins,bins);
for i=1:size(Ball,1)
    for j=1:size(Ball,2)
        H(bin1(i,j),bin2(i,j),bin3(i,j)) = H(bin1(i,j),bin2(i,j),bin3(i,j)) + 1;
    end
end

%% Read new img
img = imread('ball02.JPG');

%% Backprojecting on a new img:
new_img = BackProjection(img,H,bins,3);
%figure; imshow(new_img);

%% Filter
sigma = 1;
h = fspecial('gaussian', 5, sigma);
filtered_img = imfilter(new_img, h);
figure; imshow(filtered_img);

%% Circle finding

edged_img = edge(filtered_img,'canny');
figure; imshow(edged_img);

radiusRange = [20 40];
%[centers,radii] = imfindcircles(filtered_img,radiusRange);
[centers,radii] = imfindcircles(edged_img,radiusRange, 'Sensitivity',0.9);

%% Circle plotting
figure; imshow(img);
hold on;
viscircles(centers, radii);