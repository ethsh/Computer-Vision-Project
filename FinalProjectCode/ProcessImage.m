function ProcessImage(processing_image, originalImage)

global Config;

layer1 = double(processing_image(:,:,1));
layer1(abs(layer1)<16) = 0;

layer2 = double(processing_image(:,:,2));
layer2(abs(layer2)<10) = 0;


total_layer = Config.lc_lambda * layer1 + (1-Config.lc_lambda) * layer2;
total_layer = medfilt2(total_layer, [5 5]);
H=fspecial('gaussian',3,0.2);
total_layer=conv2(total_layer,H,'same');
p_nms = prctile(total_layer(:),Config.nms_Percentage);
% TODO : find other NMS
I_nms = NMS(total_layer,Config.minimal_distance_between_maxima,p_nms);
figure; imshow(I_nms);

processing_mat_1 = GeneratePrecessingMats(layer1, 1);
processing_mat_2 = GeneratePrecessingMats(layer2, 1);

[grad_x_1, grad_y_1] = CalcGradient(processing_mat_1, Config.grad_type);
[grad_x_2, grad_y_2] = CalcGradient(processing_mat_2, Config.grad_type);

processing_mat = Config.lc_lambda * processing_mat_1 + (1-Config.lc_lambda) * processing_mat_2;

[y, x] =ind2sub(size(I_nms), find(I_nms));
centers = [y, x, I_nms(logical(I_nms))];
centers = sortrows(centers, 3);
% centers = centers(1:floor(length(centers)),:,:);

figure; imshow(originalImage);
hold on; points_counter = 1; flag = 0;
for i = 1:floor(size(centers,1))
    y = centers(i,1);
    x = centers(i,2);
    r = Config.initial_r_guess;
    
    if (flag)
        Ezer_mat = repmat([x y],size(Mat_outputs,1),1);
        Distances_mat = Ezer_mat - Mat_outputs(:,1:2);
        distances_vec = sqrt(Distances_mat(:,1).^2 + Distances_mat(:,2).^2);
        if (sum((distances_vec<Mat_outputs(:,3))) > 0)
            % This point resides in one of the previously found circles!
            continue;
        end
    end
    
    % First Convergence - Big picture
    [x_final, y_final, r_final, S_final] = RecursiveParametersConvergence(y, x, r, processing_mat, grad_x_1, grad_y_1, grad_x_2, grad_y_2);
    if (S_final == -31000)
        points_counter = points_counter + 1;
        continue;
    end
    if (r_final <= 0)
        points_counter = points_counter + 1;
        continue;
    end
    
    % Fine tuning
    SideLength = min([round(r_final)*4+1 size(originalImage,1) size(originalImage,2)]);
    dT = (SideLength -1)/2;
    
    y_final = round(y_final);
    x_final = round(x_final);
    
    %validation
    if ((x_final - dT <= 0) || (y_final - dT <= 0) || (x_final + dT > size(layer1, 2)) || (y_final + dT > size(layer1, 1)))
        dT = min(min(x_final - 1, y_final - 1),min(size(layer1, 2) - x_final - 1,size(layer1, 1) - y_final - 1));
    end
    
    layer1_cropped = layer1(y_final-dT:y_final+dT, x_final-dT:x_final+dT);
    layer2_cropped = layer2(y_final-dT:y_final+dT, x_final-dT:x_final+dT);
    
    processing_mat_1_fine = GeneratePrecessingMats(layer1_cropped, 0);
    processing_mat_2_fine = GeneratePrecessingMats(layer2_cropped, 0);
    
    [grad_x_1_fine, grad_y_1_fine] = CalcGradient(processing_mat_1_fine, Config.grad_type);
    [grad_x_2_fine, grad_y_2_fine] = CalcGradient(processing_mat_2_fine, Config.grad_type);
    
    processing_mat_fine = Config.lc_lambda * processing_mat_1_fine + (1-Config.lc_lambda) * processing_mat_2_fine;
    if (all(all(processing_mat_fine == 0)))
        points_counter = points_counter + 1;
        continue;
    end
    
    [x_final_fine, y_final_fine, r_final_fine, S_final_fine] = RecursiveParametersConvergence(dT, dT, r_final, processing_mat_fine, grad_x_1_fine, grad_y_1_fine, grad_x_2_fine, grad_y_2_fine);
    if (S_final_fine == -31000)
        points_counter = points_counter + 1;
        continue;
    end
    
    Mat_outputs(points_counter,:) = [x_final+x_final_fine-dT, y_final+y_final_fine-dT, r_final_fine, S_final_fine];
    if (r_final>0)
        viscircles([Mat_outputs(points_counter,1),Mat_outputs(points_counter,2)], r_final_fine);
    end
    if (points_counter == 1)
        flag = 1;
    end
    points_counter = points_counter + 1;
end

end