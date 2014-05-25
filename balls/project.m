function project ()
    datestr(now)
    BB = 'MVC-00';
    for AA=1:25
        disp(strcat('Image number ',num2str(AA)));
        if (AA > 9)
            BB='MVC-0';
        end
        text=strcat(BB,num2str(AA),'F.JPG');
        pic = DoIt(text);
    end
    datestr(now)
end

function pic = DoIt (text)

    % Read Image
    warning off
    rgb = imread(strcat('../balls/',text));
   
    % Downsize to save time
    rgb_small=imresize(rgb,0.25,'bilinear');
    
    % Get RGB components seperated
    r=rgb_small(:,:,1);
    g=rgb_small(:,:,2);
    b=rgb_small(:,:,3); 
      
    % Get size of the matrixes
    [m,n]=size(r);
       
    r_bin=logical(zeros(m,n));
    g_bin=logical(zeros(m,n));
    y_bin=logical(zeros(m,n));
    
    % Color assortment - green,yellow and red
    for i=1:m
        for j=1:n
            if (((r(i,j) > 1.7*b(i,j)) && (r(i,j) > 1.7*g(i,j) && (g(i,j) > b(i,j)))) || ((r(i,j) > 1.8*g(i,j)) && (r(i,j) > 1.3*b(i,j)) && (b(i,j) > g(i,j))))
                r_bin(i,j)=1;
            end
            if (((g(i,j) > r(i,j)) && (r(i,j) > 1.4*b(i,j))) || ((g(i,j) > 1.2*r(i,j)) && (g(i,j) > 1.2*b(i,j)) && (r(i,j) > b(i,j))))
                g_bin(i,j)=1;
            end
            if ((r(i,j) > g(i,j)) && (g(i,j) > 1.4*b(i,j)) && (r(i,j) < 2*g(i,j)))
                y_bin(i,j)=1;
            end
        end
    end  
       
    % Create initialized mask matrix
    mask=logical(zeros(m,n));
    
    % Define relevant mask distance
    rad = m/12;
    
    % Check if an index (i,j) has in its neighbor all Yellow, Red & Green
    % colors as appeared on the ball
    for i=1:m
        for j=1:n
            r_mask=0;
            y_mask=0;
            g_mask=0;
            for k=max(1,i-rad):min(i+rad,m)
                for l=max(1,j-rad):min(j+rad,n)
                    if ((~r_mask)&&(r_bin(k,l) == 1))
                        r_mask=1;
                    end
                    if ((~y_mask)&&(y_bin(k,l) == 1))
                        y_mask=1;
                    end
                    if ((~g_mask)&&(g_bin(k,l) == 1))
                        g_mask=1;
                    end     
                    if (r_mask&&y_mask&&g_mask)
                        mask(i,j)=1;
                        break;
                    end
                end
                if (mask(i,j) == 1)
                    break;
                end
            end
        end
    end
    
% Combine all the components and apply the relevance mask on it
    all_bin=logical(r_bin+g_bin+y_bin);
    pic=all_bin.*mask;           
            
    % Apply median filter to fill in the remaining shapes
    pic=medfilt2(pic,[5 5]);
    
    % Remove smaller shapes than estimation of minimum ball size
    [B,L,N] = bwboundaries(pic,'noholes');
    min_perimeter = 15;
    for i=1:N
        if (i == 1)
            B_vect = length(B{i});
        else
            B_vect = [B_vect length(B{i})];
        end
    end
    if ((N > 0) && (min(B_vect) < min_perimeter))
        for i=1:m
            for j=1:n
                if ((L(i,j) > 0) && (B_vect(L(i,j)) < min_perimeter))
                    mask(i,j) = 0;
                end
            end
        end             
        % Re-calculate pic
        pic = all_bin.*mask;
        pic = medfilt2(pic,[5 5]);        
    end
        
    max_potential_circles = 4;
    treshold = 0.97;
  
    matrix_top = []; mat_1 = []; mat_2 = []; mat_3 = []; mat_4 = []; mat_5 = [];

    % Find potential circles from each size ctegory
    [centers,radii,strength] = imfindcircles(pic,[4  10],'ObjectPolarity','bright','Sensitivity',treshold,'Method','TwoStage');
    mat_1 = [strength centers radii];            
    if (numel(mat_1)/4 > max_potential_circles - 1)
        mat_1 = mat_1(1:max_potential_circles - 1,:);
    end
    [centers,radii,strength] = imfindcircles(pic,[10 20],'ObjectPolarity','bright','Sensitivity',treshold,'Method','TwoStage');
    mat_2 = [strength centers radii];
    if (numel(mat_2)/4 > max_potential_circles - 1)
        mat_2 = mat_2(1:max_potential_circles - 1,:);
    end
    [centers,radii,strength] = imfindcircles(pic,[20 30],'ObjectPolarity','bright','Sensitivity',treshold,'Method','TwoStage');
    mat_3 = [strength centers radii];
    if (numel(mat_2)/4 > max_potential_circles - 1)
        mat_3 = mat_3(1:max_potential_circles - 1,:);
    end            
    [centers,radii,strength] = imfindcircles(pic,[30 45],'ObjectPolarity','bright','Sensitivity',treshold,'Method','TwoStage');
        mat_4 = [strength centers radii];
    if (numel(mat_4)/4 > max_potential_circles - 2)
        mat_4 = mat_4(1:max_potential_circles - 2,:);
    end            
    [centers,radii,strength] = imfindcircles(pic,[45 60],'ObjectPolarity','bright','Sensitivity',treshold,'Method','TwoStage');
    mat_5 = [strength centers radii];
    if (numel(mat_5)/4 > max_potential_circles - 3)
        mat_5 = mat_5(1:max_potential_circles - 3,:);
    end      

    matrix_top = [mat_1;mat_2;mat_3;mat_4;mat_5];

    %%%%%% Leave only the most likely circles     
    if (~isempty(matrix_top))
        matrix_top = sortrows(matrix_top,-1);  
        if (length(matrix_top(:,1)) > max_potential_circles)       
            matrix_top = matrix_top(1:max_potential_circles,:);
        end

              %%% First step - Circle strength
        index=[];
        for i=1:length(matrix_top(:,1))
                if (matrix_top(1,1)/2.1) > matrix_top(i,1)
                    index(i)=i;
                end
        end
        matrix_top = matrix_top(~ismember(1:size(matrix_top, 1), index), :);
        
              %%% Second step - Missing Colors
             
        index=[];      
        for i=1:length(matrix_top(:,1))

                % yellow missing
                cy = Color_intensity(y_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))==0;
                % red missing
                cr = Color_intensity(r_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))==0;
                % green missing
                cg = Color_intensity(r_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))==0;

                if (length(matrix_top(:,1))>1)&&(cy||cg||cr)
                index(i)=i;
                end
        end
        matrix_top = matrix_top(~ismember(1:size(matrix_top, 1), index), :);
        
               %%% Third step - Color combinations
               
        index=[];      
        for i=1:length(matrix_top(:,1))

                % too much yellow
                cy_high = Color_intensity(y_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))>((pi*matrix_top(i,4)^2)/2.3);
                % Huge amount of yellow
                cy_huge = Color_intensity(y_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))>((pi*matrix_top(i,4)^2)/1.6);
                % too much red
                cr_high = Color_intensity(r_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))>((pi*matrix_top(i,4)^2)/2);
                % not enough green
                cg_low = Color_intensity(g_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))<((pi*matrix_top(i,4)^2)/5);
                % not enough red
                cr_low = Color_intensity(r_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))<((pi*matrix_top(i,4)^2)/8);
                % not enough green
                cg_low2 = Color_intensity(g_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))<((pi*matrix_top(i,4)^2)/8);
                % not enough yellow
                cy_low = Color_intensity(y_bin, round(matrix_top(i,2)), round(matrix_top(i,3)), matrix_top(i,4))<((pi*matrix_top(i,4)^2)/8);
                
                if (length(matrix_top(:,1))>1)&&((cy_high&&cr_low&&cg_low2)||(cy_high&&cg_low2)||(cr_high&&cg_low2&&cy_low)||(cy_huge&&cg_low))
                index(i)=i;
                end
        end
        matrix_top = matrix_top(~ismember(1:size(matrix_top, 1), index), :);
        
        %%% Final step - Merge simmilar objects
    
        if length(matrix_top(:,1))>1
            duplicates = true;
            while duplicates
                duplicates = false;
                for i=1:length(matrix_top(:,1))-1
                    for k=i+1:length(matrix_top(:,1))
                        circ_dist = pdist([matrix_top(i,2),matrix_top(i,3);matrix_top(k,2),matrix_top(k,3)],'euclidean');
                        if abs(circ_dist)<(matrix_top(i,4)+matrix_top(k,4))
                           if  matrix_top(i,4)>matrix_top(k,4)
                               matrix_top = matrix_top(~ismember(1:size(matrix_top, 1), k), :);
                           else
                               matrix_top(i,:) = matrix_top(k,:);
                               matrix_top = matrix_top(~ismember(1:size(matrix_top, 1), k), :);
                           end
                           duplicates = true;
                           break
                        end
                    end
                    if duplicates 
                        break 
                    end
                end
            end 
        end
    
    end
    
    % Load image
    handle = figure ('visible','off'); imshow(rgb);  
    
    % Show circles if such exists after transforming back to original scale
    if (~isempty(matrix_top))     
        y = viscircles((diag([4 4])*matrix_top(:,2:3)')', matrix_top(:,4)*4+2,'LineWidth',4);
    end
    
    % Save result
    print(handle, '-djpeg', strcat('../results/',text));
end

function m = Color_intensity(img, xc, yc, r)
% Color_intensity computes number of values inside a circle
% This section is for efficiency only - avoids wasting computation time on
% pixels outside the bounding square
[sy sx] = size(img);
xmin = max(1, floor(xc-r));
xmax = min(sx, ceil(xc+r));
ymin = max(1, floor(yc-r));
ymax = min(sy, ceil(yc+r));
img = img(ymin:ymax, xmin:xmax); % trim boundaries
xc = xc - xmin + 1;
yc = yc - ymin + 1;
% Make a circle mask
[x y] = meshgrid(1:size(img,2), 1:size(img,1));
mask = (x-xc).^2 + (y-yc).^2 < r.^2;
% Compute color intensity
m = sum(sum(double(img) .* mask));
end
