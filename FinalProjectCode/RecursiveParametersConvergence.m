function [x, y, r] = RecursiveParametersConvergence(y0, x0, r0, I, I1, I2, grad_x_1, grad_y_1, grad_x_2, grad_y_2, last_integral, counter)

global Config;

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


if counter == Config.recursive_max_iter
    y = y0;
    x = x0;
    r = r0;
    return;
end


mat_ezer = zeros(2*r0+1);
[x_grid, y_grid] = meshgrid(-r0:r0,-r0:r0);
mat_ezer((x_grid.^2 + y_grid.^2) <= r0^2) = 1;
mat_ezer_r = sqrt(x_grid.^2 + y_grid.^2);
mat_ezer_r((x_grid.^2 + y_grid.^2) > r0^2) = 0;
relevant_I = I(y0-r0:y0+r0, x0-r0:x0+r0);
S = sum(sum(double(relevant_I).*mat_ezer));

if (last_integral > 0)
    integral_diff = abs((S-last_integral)/last_integral);
else
    integral_diff = S/Config.initial_integral_diff_facotr;
end

if abs(integral_diff) < Config.error_allowed
    y = y0;
    x = x0;
    r = r0;
    return;
end

relevant_grad_x_1 = grad_x_1(y0-r0:y0+r0, x0-r0:x0+r0);
relevant_grad_y_1 = grad_y_1(y0-r0:y0+r0, x0-r0:x0+r0);

relevant_grad_x_2 = grad_x_2(y0-r0:y0+r0, x0-r0:x0+r0);
relevant_grad_y_2 = grad_y_2(y0-r0:y0+r0, x0-r0:x0+r0);

relevant_grad_x = Config.lc_lambda * relevant_grad_x_1 + (1 - Config.lc_lambda) * relevant_grad_x_2;
relevant_grad_y = Config.lc_lambda * relevant_grad_y_1 + (1 - Config.lc_lambda) * relevant_grad_y_2;

[dy, dx, dr] = CalcGradSteps(relevant_I, relevant_grad_x, relevant_grad_y, mat_ezer, mat_ezer_r, r0); % sending mat_ezer so we won't have to calculate it again inside the function

steps = [dy, dx, dr];
nor = norm(steps);
% steps = -steps/norm(steps)*Config.step_size_factor;
steps = -steps/nor*3;

disp(counter)
disp(r0)
disp(x0)
disp(y0)
disp(S)

if (steps(3) == inf || isnan(steps(3)))
    y = y0;
    x = x0;
    r = r0;
    return;
end

[y, x, r] = RecursiveParametersConvergence(y0 + steps(1), x0 + steps(2), r0 + steps(3), I,I1, I2, grad_x_1, grad_y_1, grad_x_2, grad_y_2, S, counter+1);

end