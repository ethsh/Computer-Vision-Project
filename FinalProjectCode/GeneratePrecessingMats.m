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

% p_initial = prctile(output_mat(:),zero_percent);
% output_mat(output_mat<=p_initial-Config.prctile_security) = 0;

output_mat = (output_mat-mean(output_mat(:))/(std(((output_mat(:))))));

%output_mat = medfilt2(output_mat, [20 20]);

se = strel('disk',Config.dilate_len);
% output_mat = imdilate(output_mat,se);
try
output_mat = imdilate(output_mat,se);
catch
    aaa = 1;
end
% h = fspecial('gaussian', 10, 1);
% output_mat = conv2(output_mat, h, 'same');

output_mat = output_mat/max(output_mat(:));
p_initial = (min(output_mat(:)));
%p_initial = prctile(output_mat(:),penatly_percent);
output_mat(output_mat <= 0) = penalty;
%output_mat(output_mat == 0) = penalty;

end