function [dy, dx, dr] = CalcGradSteps(relevant_I, relevant_grad_x, relevant_grad_y, mat_ezer_r, mat_ezer_empty, r0)

dy = sum(sum(relevant_grad_y.*mat_ezer_r));
dx = sum(sum(relevant_grad_x.*mat_ezer_r));
dr = -sum(sum(r0*double(relevant_I).*mat_ezer_empty))/1.2;

end