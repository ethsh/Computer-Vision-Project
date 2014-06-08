function output_mat = GeneratePrecessingMats(input_mat, initial_flag)

global Config

if (initial_flag)
    zero_percent = Config.zero_Percentage;
    penatly_percent = Config.penalty_Percentage;
    penalty = Config.penalty;
else
    zero_percent = Config.fine_zero_Percentage;
    penatly_percent = Config.fine_penalty_Percentage;
    penalty = Config.fine_penalty;
end

output_mat = abs(input_mat);
%output_mat(output_mat < Config.initial_threshold) = 0;
p_initial = prctile(output_mat(:),zero_percent);
output_mat(output_mat<=p_initial) = 0;

output_mat = (output_mat-mean(output_mat(:))/(std(std((output_mat)))));

%output_mat = medfilt2(output_mat, [20 20]);

se = strel('disk',Config.dilate_len);
% output_mat = imdilate(output_mat,se);
output_mat = imdilate(output_mat,se);

% h = fspecial('gaussian', 10, 1);
% output_mat = conv2(output_mat, h, 'same');

output_mat = output_mat/max(output_mat(:));
p_initial = prctile(output_mat(:),penatly_percent);
output_mat(output_mat <= p_initial) = penalty;

% figure; surf(output_mat);
% figure; imshow(output_mat);
return;


output_mat=abs(input_mat);
T=prctile(output_mat(1:end),zero_percent);
output_mat(output_mat<T)=0;
output_mat=(output_mat-mean(output_mat(:)))./(std(output_mat(:)));

se = strel('disk',Config.dilate_len);
output_mat = imdilate(output_mat,se);

output_mat=(output_mat./max(output_mat(:)));

T=prctile(output_mat(1:end),penatly_percent);
output_mat(output_mat<=T)=penalty;

h = fspecial('gaussian', 10, 1);
output_mat = conv2(output_mat, h, 'same');


% figure; surf(output_mat);
% figure; imshow(output_mat);

end