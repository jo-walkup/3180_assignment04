function [error,h_avg, num_evals] = global_error(rate_func_in, sol_list, XA, t_range, href, BT)

% function [error,h_avg,num_evals] = global_error(rate_func_in,sol_func,XA,t_range,href,BT)
% 
%     [~,X_list,h_avg,num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_range, XA, href, BT);
%     X_numerical_final = X_list(end,:)';
%     X_exact_final = sol_func(t_range(end));
%     error = norm(X_numerical_final - X_exact_final);
% end
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,t_range,XA,href,BT);; % change explicit_midpoint_step to the other funcs when needed
    error = norm(X_list(end,:)' - sol_list(t_range(end)));
end