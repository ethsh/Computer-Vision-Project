function Main()

close all;
clear parameters;
clc;

%% configuration
HSV_flag = 0; % 1 for HSV, 0 for Lab
if HSV_flag
    wanted_dimensions = 1:2;
else
    wanted_dimensions = 2:3;
end

global Config;

Config.initial_threshold = 16;
Config.zero_Percentage = 90;
Config.penalty_Percentage = 95;
Config.fine_zero_Percentage = 50;
Config.fine_penalty_Percentage = 50;
Config.fine_penalty= -1;


Config.nms_Percentage = 99;
Config.penalty= -1;
Config.dilate_len = 5;

Config.lc_lambda = 0.5;


Config.minimal_distance_between_maxima = 5;
Config.grad_type = 'prewitt';
Config.initial_r_guess = 65;
Config.recursive_max_iter = 30;
Config.error_allowed = 10^-2;
%Config.step_size_factor = 0.005;
%Config.r_additional_step_size_facotr = 0.05;
%Config.initial_integral_diff_facotr = 100;

disp('Begining...');

for i=1:25
    close all;
    disp(['Image ' num2str(i)]);
    tic;
    %% Read new img
    if (i<10)
        num = ['0' num2str(i)];
    else
        num = num2str(i);
    end
    file = ['balls/MVC-0' num 'F.JPG'];
    [img,map_img] = imread(file);
    figure; imshow(img);
    
    wanted_image = extract_wanted_format(img, HSV_flag);
    image_working_dimensions = wanted_image(:,:,wanted_dimensions);
    
    ProcessImage(image_working_dimensions, img);
    
    toc;
end

end

function output_img = extract_wanted_format(rgb_img, hsv_flag)
if (hsv_flag)
    output_img = rgb2hsv(rgb_img);
else
    [output_img(:,:,1), output_img(:,:,2), output_img(:,:,3)] = RGB2Lab(rgb_img(:,:,1),rgb_img(:,:,2),rgb_img(:,:,3));
%     
%     colorTransform = makecform('srgb2lab');
%     output_img = applycform(rgb_img, colorTransform);
end
end