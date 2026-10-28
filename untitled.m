%example for how to use compute_planetary_motion(...)
function usage_example_comparisons()
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    clf();
    clear all;

   % Dormand-Prince Butcher Tableau

DormandPrince = struct();

DormandPrince.C = ...
    [0, 1/5, 3/10, 4/5, 8/9, 1, 1];

DormandPrince.B = ...
    [35/384, 0, 500/1113, 125/192, ...
     -2187/6784, 11/84, 0;
     5179/57600, 0, 7571/16695, 393/640, ...
     -92097/339200, 187/2100, 1/40];

DormandPrince.A = ...
    [0, 0, 0, 0, 0, 0, 0;
     1/5, 0, 0, 0, 0, 0, 0;
     3/40, 9/40, 0, 0, 0, 0, 0;
     44/45, -56/15, 32/9, 0, 0, 0, 0;
     19372/6561, -25360/2187, 64448/6561, ...
        -212/729, 0, 0, 0;
     9017/3168, -355/33, 46732/5247, ...
        49/176, -5103/18656, 0, 0;
     35/384, 0, 500/1113, 125/192, ...
        -2187/6784, 11/84, 0];

    
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 5;
    y0 = 3;
    dxdt0 = 0;
    dydt0 = 1.8;
    

    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    h_ref=0.0005;

    G_rate= @(t,v) gravity_rate_func(t, v,orbit_params);
    p=5.5;
    error_desired=0.0001;
    
    hs = logspace(-3,-1,100);

    num_evals_all = zeros(size(hs));
    failed_steps_all = zeros(size(hs));

    for n = 1:length(hs)
        h = hs(n);
        [t_list,X_list,h_avg, num_evals failed_steps] = explicit_RK_variable_step_integration(G_rate,t_range,V0,h_ref,DormandPrince,p,error_desired);
        num_evals_all(n) = num_evals;
        failed_steps_all(n) = failed_step;

    end

    figure();
    loglog(h_avg_list_fe, global_errors_fe, "b.", MarkerSize=20); hold on
    coeffs_fe = polyfit(log(hs), log(global_errors_fe), 1);
    globalfit_fe = exp(polyval(coeffs_fe, log(hs)));
    loglog(hs, globalfit_fe, "r", LineWidth=3);

    loglog(h_avg_list_em, global_errors_em, "m.", MarkerSize=10); hold on
    coeffs_em = polyfit(log(hs), log(global_errors_em), 1);
    globalfit_em = exp(polyval(coeffs_em, log(hs)));
    loglog(hs, globalfit_em, "g", LineWidth=1);

    loglog(h_avg_list_rk, global_errors_rk, "k.", MarkerSize=5); hold on
    coeffs_rk = polyfit(log(hs), log(global_errors_rk), 1);
    globalfit_rk = exp(polyval(coeffs_rk, log(hs)));
    loglog(hs, globalfit_rk, "c", LineWidth=1);
    
    title("Global Truncation Error for Various Runge-Kutta Methods",'FontSize',13);
    xlabel("$h_{avg} (-)$ ",'FontSize',13);
    ylabel("Global Truncation error $(-)$",'FontSize',13)
    legend("Forward Euler error", "Forward Euler fit", "Implicit midpoint error", "Implicit midpoint fit", "RK4 error", "RK4 fit", 'FontSize',11, "Location", "northwest")

    myfig = gcf();
    exportgraphics(myfig , 'global_err_explicit_01.png', 'Resolution', 600);
end
end