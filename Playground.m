close all;
clear all;
clc;

sigma = 2;
radiusRange = [40 80];

[im, map] = imread('ball02.jpg');
r = 174;
g = 121;
b = 139;
diff = 4;

a = zeros(size(im));

find(im(:,:,1)>170 && im(:,:,1)<178 && im(:,:,2)>117 && im(:,:,2)<125 && im(:,:,3)>135 && im(:,:,3)<143);


figure; imshow(im);
im = rgb2hsv(im);
imshow(im)

% h = fspecial('gaussian', 20, sigma);
% im = imfilter(im, h);
% figure; imshow(im);

dr = edge(im(:,:,1),'canny', [], sigma);
dg = edge(im(:,:,2),'canny', [], sigma);
db = edge(im(:,:,3),'canny', [], sigma);

figure; imshow(dr); % 1
figure; imshow(dg); % 1
figure; imshow(db); % 1

% dr = edge(im(:,:,1));
% dg = edge(im(:,:,2));
% db = edge(im(:,:,3));
% 
% figure; imshow(dr); % 1
% figure; imshow(dg); % 1
% figure; imshow(db); % 1

d = double(dr) + double(dg) + double(db);
a = zeros(size(d));
a(d>=2) = 1;

figure; imshow(a); % 1

[centers,radii] = imfindcircles(100*a,radiusRange, 'Sensitivity', 0.98);

imshow(im);
hold on;
viscircles(centers, radii)



[centersR,radiiR] = imfindcircles(dr,radiusRange);
[centersG,radiiG] = imfindcircles(dg,radiusRange);
[centersB,radiiB] = imfindcircles(db,radiusRange);

[centersTot, radiiTot] = imfindcircles(im2,radiusRange,'Sensitivity',0.9)

theta = 0:0.05:2*pi;
for j=1:size(centersR, 1)
    x = centersR(j, 1)+radiiR*sin(theta);
    y = centersR(j, 2)+radiiR*cos(theta);
    for i=1:length(x)
        im(round(y(i)),round(x(i)),:) = 255;
    end
end

for j=1:size(centersG, 1)
    x = centersG(j, 1)+radiiG*sin(theta);
    y = centersG(j, 2)+radiiG*cos(theta);
    for i=1:length(x)
        im(round(y(i)),round(x(i)),:) = 255;
    end
end

for j=1:size(centersB, 1)
    x = centersB(j, 1)+radiiB*sin(theta);
    y = centersB(j, 2)+radiiB*cos(theta);
    for i=1:length(x)
        im(round(y(i)),round(x(i)),:) = 255;
    end
end



figure; imshow(im);
%%

im2 = conv2(h, im);
figure; imshow(im, map); % 1

BW = edge(im,'canny', [], 10);

im2 = im2(10:end-10, 10:end-10);
BW = BW(10:end-10, 10:end-10);

figure; imshow(BW); % 2

thres = nan(size(im2));
thres(BW) = im2(BW);
thres = inpaint_nans(thres);

output = 255*ones(size(im2));
output(im>thres/2.2) = 0;

figure; imshow(output); %3

output = conv2(h, output);
output = output/max(output(:));
output = round(output);
output = output*255;

figure; imshow(output); %4