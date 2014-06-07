close all; clear all; clc;

[img,map_img] = imread('balls/MVC-005F.JPG');
std_r = 20;
std_l = std_r/2;
w = size(img, 1);
h = size(img, 2);
r_std = zeros(w,h);
g_std = zeros(w,h);
b_std = zeros(w,h);

square = std_l*[1 1];

for i=std_l+1:w-std_l-1
    for j=std_l+1:h-std_l-1
        temp = impart(img,[j, i],square);
        r_std(i, j) = std(reshape(double(temp(:, :, 1)),1,[]));
        g_std(i, j) = std(reshape(double(temp(:, :, 2)),1,[]));
        b_std(i, j) = std(reshape(double(temp(:, :, 3)),1,[]));
    end
end

figure; imshow(r_std, map_img);
figure; imshow(g_std, map_img);
figure; imshow(b_std, map_img);