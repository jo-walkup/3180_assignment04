function [error,h_avg, num_evals] = global_error(rate_func_in, sol_func, XA, t_range, href, BT)
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(G_rate,t_range,V0,h_ref,BT_fe);; % change explicit_midpoint_step to the other funcs when needed
    error = norm(X_list(:,end) - sol_func(t_range(end)));
end