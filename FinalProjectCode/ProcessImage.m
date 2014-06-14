function ProcessImage(processing_image, originalImage)

global Config;

figure; imshow(originalImage);

original_processing_image = processing_image;
processing_image = imresize(processing_image, Config.scale_factor);

layer1 = double(processing_image(:,:,1));
layer1(abs(layer1)<Config.initial_threshold) = 0;

layer2 = double(processing_image(:,:,2));
layer2(abs(layer2)<Config.initial_threshold) = 0;


total_layer = Config.lc_lambda * layer1 + (1-Config.lc_lambda) .* layer2;
total_layer = medfilt2(total_layer, [5 5]);
H=fspecial('gaussian',3,0.2);
total_layer=conv2(total_layer,H,'same');
p_nms = prctile(total_layer(:),Config.nms_Percentage);
% TODO : find other NMS
I_nms = NMS(total_layer,Config.minimal_distance_between_maxima,p_nms);
%figure; imshow(I_nms);

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
    
    % convert to original values
    
    x_final = round(x_final/Config.scale_factor);
    y_final = round(y_final/Config.scale_factor);
    r_final = round(r_final/Config.scale_factor);
    
    layer1 = double(original_processing_image(:,:,1));
    layer1(abs(layer1)<16) = 0;
    
    layer2 = double(original_processing_image(:,:,2));
    layer2(abs(layer2)<10) = 0;
    
    % Fine tuning
    SideLength = min([round(r_final)*4+1 size(originalImage,1) size(originalImage,2)]);
    dT = (SideLength -1)/2;
    
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

    if (points_counter == 1)
        flag = 1;
    end
    points_counter = points_counter + 1;
end

Mat_outputs(find(Mat_outputs(:,1) == 0),:)=[];
Mat_outputs_final = [Mat_outputs(1,:) 1]; k = 1;
for i=2:size(Mat_outputs,1)
    bool = true;
    for j=1:size(Mat_outputs_final,1)
        if (sum(abs(Mat_outputs(i,1:3)-Mat_outputs_final(j,1:3)) < [10 10 10]) == 3)
            bool = false;
            break;
        end
    end
    if (bool)
        k=k+1;
        Mat_outputs_final(k,1:4) = Mat_outputs(i,1:4);
        Mat_outputs_final(k,5) = 1;
    else
        Mat_outputs_final(j,1:4) = (Mat_outputs_final(j,5)*Mat_outputs_final(j,1:4) + Mat_outputs(i,1:4))/(Mat_outputs_final(j,5)+1);
        Mat_outputs_final(j,5) = Mat_outputs_final(j,5)+1;
    end        
end

figure; imshow(originalImage);
hold on;
for i=1:size(Mat_outputs_final,1)
    
    x_final = round(Mat_outputs_final(i,1));
    y_final = round(Mat_outputs_final(i,2));
    r_final = round(Mat_outputs_final(i,3));
    
    SideLength = min([round(r_final*3) size(originalImage,1) size(originalImage,2)]);
    dT = (SideLength -1)/2;
    
    %validation
    if ((x_final - dT <= 0) || (y_final - dT <= 0) || (x_final + dT > size(layer1, 2)) || (y_final + dT > size(layer1, 1)))
        dT = min(min(x_final - 1, y_final - 1),min(size(layer1, 2) - x_final - 1,size(layer1, 1) - y_final - 1));
    end
    
    image_cropped = originalImage(y_final-dT:y_final+dT, x_final-dT:x_final+dT, :);
    image_cropped = rgb2gray(image_cropped);
    
%     H=fspecial('gaussian',5,1);
%     image_cropped=conv2(double(image_cropped),H,'same');
    I = edge(image_cropped,'canny');
    [centers,radii] = imfindcircles(I,round([0.8*r_final 1.2*r_final]),'Sensitivity',0.96);
    
    if (~isempty(centers))
        for m=1:size(centers,1)
            x_descale = x_final - dT + centers(m,1);
            y_descale = y_final - dT + centers(m,2);
            r_descale = radii(m);
            
            viscircles([x_descale,y_descale], r_descale);
        end
    end
end



% for i=1:size(Mat_outputs_final,1)
%     if (Mat_outputs_final(i,3)>0)
%         viscircles([Mat_outputs_final(i,1),Mat_outputs_final(i,2)], Mat_outputs_final(i,3));
%     end
% end

end