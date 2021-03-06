function [S, steps] = ParametersEstimation(y0, x0, r0, I, grad_x_1, grad_y_1, grad_x_2, grad_y_2)

if (r0==0)
    S = 0;
    steps_normed = [0;0;0];
    return;
end

global Config;

% Inputs validation
y0 = round(y0);
x0 = round(x0);
r0 = round(r0);

if (x0 - r0 <= 0)
    r0 = x0 - 1;
end
if (y0 - r0 <= 0)
    r0 = y0 - 1;
end
if (x0 + r0 > size(I, 2))
    r0 = size(I, 2) - x0 - 1;
end
if (y0 + r0 > size(I, 1))
    r0 = size(I, 1) - y0 - 1;
end

% Matrices in the relevant area of the circle:
mat_ezer = zeros(2*r0+1);
[x_grid, y_grid] = meshgrid(-r0:r0,-r0:r0);
mat_ezer((x_grid.^2 + y_grid.^2) <= r0^2) = 1;

mat_ezer_temp = zeros(2*r0+1);
mat_ezer_temp((x_grid.^2 + y_grid.^2) <= (r0-1)^2) = 1;

mat_ezer_empty = mat_ezer-mat_ezer_temp;

mat_ezer_r = sqrt(x_grid.^2 + y_grid.^2);
mat_ezer_r((x_grid.^2 + y_grid.^2) > r0^2) = 0;
try
relevant_I = I(y0-r0:y0+r0, x0-r0:x0+r0);
catch
    aaa =1;
end
% Integral + Advance Calculation                                           
S = -sum(sum(double(relevant_I).*mat_ezer_r));

relevant_grad_x_1 = grad_x_1(y0-r0:y0+r0, x0-r0:x0+r0);
relevant_grad_y_1 = grad_y_1(y0-r0:y0+r0, x0-r0:x0+r0);

relevant_grad_x_2 = grad_x_2(y0-r0:y0+r0, x0-r0:x0+r0);
relevant_grad_y_2 = grad_y_2(y0-r0:y0+r0, x0-r0:x0+r0);

relevant_grad_x = Config.lc_lambda * relevant_grad_x_1 + (1 - Config.lc_lambda) * relevant_grad_x_2;
relevant_grad_y = Config.lc_lambda * relevant_grad_y_1 + (1 - Config.lc_lambda) * relevant_grad_y_2;

[dy, dx, dr] = CalcGradSteps(relevant_I, relevant_grad_x, relevant_grad_y, mat_ezer_r, mat_ezer_empty, r0); % sending mat_ezer so we won't have to calculate it again inside the function

steps = -[dy, dx, dr];

if (norm(steps) == 0)
    y = y0;
    x = x0;
    r = r0;
    return;
end

if (abs(steps(3)) == inf || isnan(steps(3)))
    y = y0;
    x = x0;
    r = r0;
    return;
end

end