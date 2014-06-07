function ProcessImage(processing_image, originalImage)

global Config;

layer1 = double(processing_image(:,:,1));
layer2 = double(processing_image(:,:,2));

total_layer = Config.lc_lambda * layer1 - (1-Config.lc_lambda) * layer2;
total_layer = medfilt2(total_layer, [20 20]);
H=fspecial('gaussian',3,1);
total_layer=conv2(total_layer,H,'same');
p_nms = prctile(total_layer(:),Config.nms_Percentage);
% TODO : find other NMS
I_nms = NMS(total_layer,Config.minimal_distance_between_maxima,p_nms);
figure; imshow(I_nms);

processing_mat_1 = GenerateCostMap( layer1,11,Config.penalty,Config.zero_Percentage,Config.penalty_Percentage,1);
% processing_mat_1 = GeneratePrecessingMats(layer1, 1);
processing_mat_2 = GenerateCostMap( layer2,11,Config.penalty,Config.zero_Percentage,Config.penalty_Percentage,1);
% processing_mat_2 = GeneratePrecessingMats(layer2, 1);



[grad_x_1, grad_y_1] = CalcGradient(processing_mat_1, Config.grad_type);
[grad_x_2, grad_y_2] = CalcGradient(processing_mat_2, Config.grad_type);


processing_mat = Config.lc_lambda * processing_mat_1 + (1-Config.lc_lambda) * processing_mat_2;

[y, x] =ind2sub(size(I_nms), find(I_nms));
centers = [y, x, I_nms(logical(I_nms))];
centers = sortrows(centers, 3);
% centers = centers(1:floor(length(centers)),:,:);

figure; imshow(originalImage);
hold on;
for i = 1:size(centers,1)
    y = centers(i,1);
    x = centers(i,2);
    r = Config.initial_r_guess;
    
    [x_final, y_final, r_final] = RecursiveParametersConvergence(y, x, r, processing_mat,processing_mat_1, processing_mat_2, grad_x_1, grad_y_1, grad_x_2, grad_y_2, -1, 0);
    
    
%     processing_mat_1 = GeneratePrecessingMats(double(processing_image(:,:,1)), 0);
%     processing_mat_2 = GeneratePrecessingMats(double(processing_image(:,:,2)), 0);
%     [grad_x_1, grad_y_1] = CalcGradient(processing_mat_1, Config.grad_type);
%     [grad_x_2, grad_y_2] = CalcGradient(processing_mat_2, Config.grad_type);
%     
%     
%     processing_mat = Config.lc_lambda * processing_mat_1 + (1-Config.lc_lambda) * processing_mat_2;
% 
%     
%     [x_final, y_final, r_final] = RecursiveParametersConvergence(y_final, x_final, r_final, processing_mat,processing_mat_1, processing_mat_2, grad_x_1, grad_y_1, grad_x_2, grad_y_2, -1, 0);
    if (r_final>0)
        viscircles([x_final, y_final], r_final);
    end
end
close all;
end