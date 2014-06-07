function [grad_x, grad_y] = CalcGradient(I, grad_type)

if strcmp(grad_type, 'sobel')
    g_x = [-1 0 1; -2 0 2; -1 0 1];
    g_y = g_x';
elseif strcmp(grad_type, 'prewitt')
    g_x = 1/6 * [-1 0 1; -1 0 1; -1 0 1];
    g_y = g_x';
end

grad_x = conv2(I, g_x, 'same');
grad_y = conv2(I, g_y, 'same');
end