function [x, y, r, S] = RecursiveParametersConvergence(y0, x0, r0, I, grad_x_1, grad_y_1, grad_x_2, grad_y_2)
global Config;

initial_step_size = 1* r0 / 7;
[S_old, steps] = ParametersEstimation(y0,x0,r0,I,grad_x_1,grad_y_1,grad_x_2,grad_y_2);

counter = 0; partial_change = inf; % initialization.
while ((counter < Config.recursive_max_iter) && (partial_change > Config.error_allowed))
    y = y0 + initial_step_size * steps(1)/norm(steps);
    x = x0 + initial_step_size * steps(2)/norm(steps);
    r = r0 + initial_step_size * steps(3)/norm(steps);
    
    if ((r<=0) || r == inf || isnan(r))
        y = y0;
        x = x0;
        r = r0;
        S = -31000;
        return;
    end
    
    [S_new, steps] = ParametersEstimation(y,x,r,I,grad_x_1,grad_y_1,grad_x_2,grad_y_2);
    if (norm(steps) == 0)
        y = y0;
        x = x0;
        r = r0;
        S = -31000;
        return;
    end
    
    step = initial_step_size;
    steps_old = steps;
    while (S_new > S_old + 0.05*step*norm(steps) && step>1)
        step = 0.5*step;
        y = y0 + step * steps_old(1)/norm(steps_old);
        x = x0 + step * steps_old(2)/norm(steps_old);
        r = r0 + step * steps_old(3)/norm(steps_old);
        
        try
        [S_new, steps] = ParametersEstimation(y,x,r,I,grad_x_1,grad_y_1,grad_x_2,grad_y_2);
        catch
            aaa=1;
        end
    end
    if (norm(steps) == 0)
        y = y0;
        x = x0;
        r = r0;
        S = -31000;
        return;
    end
    
    y0 = y; x0 = x; r0 = r;
    counter = counter + 1;
    partial_change = abs((S_new - S_old) / S_old);
    
    S_old = S_new;
end

S = S_new;

disp(r0)
disp(x0)
disp(y0)
disp(S_new)

end