function [dy, dx, dr] = CalcGradSteps(relevant_I, relevant_grad_x, relevant_grad_y, mat_ezer, mat_ezer_r, r0)

dx = sum(sum(relevant_grad_x.*mat_ezer));
dy = sum(sum(relevant_grad_y.*mat_ezer));

filt = [0 -1 0; -1 4 -1; 0 -1 0]/4;
mat_ezer_2 = conv2(mat_ezer,filt, 'same');
mat_ezer_2(logical(mat_ezer_2)) = 1;
dr = -sum(sum(r0*double(relevant_I).*mat_ezer_2))/4;

% dr = sum(sum(r0*double(relevant_I(1:end-1, 1:end)).*abs(diff(mat_ezer))));

end