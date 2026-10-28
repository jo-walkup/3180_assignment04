function error_calcs()

    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    clf();
    clear all;

    BT_fe=struct();

    BT_fe.A=[0];
    BT_fe.B=[1];
    BT_fe.C=[0];

    BT_rk=struct();

    BT_rk.A=[0, 0, 0, 0; 1/2, 0, 0, 0; 0, 1/2, 0, 0; 0, 0, 1, 0];
    BT_rk.B=[1/6; 1/3; 1/3; 1/6];
    BT_rk.C=[0; 0.5; 0.5; 1];

    BT_em=struct();

    BT_em.A=[0, 0; 0.5, 0];
    BT_em.B=[0;1];
    BT_em.C=[0; 0.5];

    hs = logspace(-3.5,-1,10);

    global_errors_fe = zeros(size(hs));
    global_errors_em = zeros(size(hs));
    global_errors_rk = zeros(size(hs));

    h_avg_list_fe = zeros(size(hs));
    h_avg_list_em = zeros(size(hs));
    h_avg_list_rk = zeros(size(hs));

    time_int = [0, 30];

    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;

    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    h_ref=0.005;

    G_rate= @(t,v) gravity_rate_func(t, v,orbit_params);
    
    for n = 1:length(hs)
        h = hs(n);
        [g_efe,h_avgfe, evalfe] = global_error(G_rate, V_list, V0, time_int, h, BT_fe);
        global_errors_fe(n) = g_efe;
        h_avg_list_fe(n) = h_avgfe;

        [g_eem,h_avgem, evalem] = global_error(G_rate, V_list, V0, time_int, h, BT_em);
        global_errors_em(n) = g_eem;
        h_avg_list_em(n) = h_avgem;
% compute_planetary_motion(time_int(1), V0, orbit_params)
        [g_erk,h_avgrk, evalrk] = global_error(G_rate, V_list, V0, time_int, h, BT_rk);
        global_errors_rk(n) = g_erk;
        h_avg_list_rk(n) = h_avgrk;
    end
    global_errors_fe
    global_errors_em
    global_errors_rk

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