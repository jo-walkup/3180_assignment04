function error_calcs()

    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    clf();
    clear all;

    % ralston
    BT_fe=struct();

    BT_fe.A=[0, 0; 2/3, 0];
    BT_fe.B=[1/4; 3/4];
    BT_fe.C=[0; 2/3];

    BT_rk=struct();

    BT_rk.A=[0, 0, 0, 0; 1/2, 0, 0, 0; 0, 1/2, 0, 0; 0, 0, 1, 0];
    BT_rk.B=[1/6; 1/3; 1/3; 1/6];
    BT_rk.C=[0; 0.5; 0.5; 1];

    BT_em=struct();

    BT_em.A=[0, 0; 0.5, 0];
    BT_em.B=[0;1];
    BT_em.C=[0; 0.5];

    BT_ny = struct();

    BT_ny.A = [0, 0, 0, 0, 0, 0;...
               1/3, 0, 0, 0, 0, 0;...
               4/25, 6/25, 0, 0, 0, 0;...
               1/4, -3, 15/4, 0, 0, 0;...
               2/27, 10/9, -50/81, 8/81, 0, 0;...
               2/25, 12/25, 2/15, 8/75, 0, 0];
    BT_ny.B = [23/192; 0; 125/192; 0; -27/64; 125/192];
    BT_ny.C = [0; 1/3; 2/5; 1; 2/3; 4/5];

    hs = logspace(-3,0,51);
    

    global_errors_fe = zeros(size(hs));
    global_errors_em = zeros(size(hs));
    global_errors_rk = zeros(size(hs));

    h_avg_list_fe = zeros(size(hs));
    h_avg_list_em = zeros(size(hs));
    h_avg_list_rk = zeros(size(hs));

    tstart = 0;
    tend = 3;

    time_int = [tstart, tend];

    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;

    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(tstart,tend,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    h_ref=0.005;

    G_rate= @(t,v) gravity_rate_func(t, v,orbit_params);
    
    for n = 1:length(hs)
        h = hs(n);
        [g_efe,h_avgfe, evalfe] = global_error(G_rate, V_list, V0, time_int, h, BT_fe);
        global_errors_fe(n) = g_efe;
        h_avg_list_fe(n) = h_avgfe;

        [g_eem,h_avgem, evalem] = global_error(G_rate, V_list, V0, time_int, h, BT_ny);
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
    coeffs_em = polyfit(log(hs(23:end)), log(global_errors_em(23:end)), 1); % THIS FILTERING IS SO CHOPPED!!!! MANUALLY PICKING THE LINEAR SECTION BUT WHO GAF
    globalfit_em = exp(polyval(coeffs_em, log(hs)));
    loglog(hs, globalfit_em, "g", LineWidth=1);

    loglog(h_avg_list_rk, global_errors_rk, "k.", MarkerSize=5); hold on
    coeffs_rk = polyfit(log(hs(13:end)), log(global_errors_rk(13:end)), 1);
    globalfit_rk = exp(polyval(coeffs_rk, log(hs)));
    loglog(hs, globalfit_rk, "c", LineWidth=1);
    
    title("Global Truncation Error for Various Runge-Kutta Methods",'FontSize',13);
    xlabel("$h_{avg} (-)$ ",'FontSize',13);
    ylabel("Global Truncation error $(-)$",'FontSize',13)
    legend("Ralston error", "Ralston fit", "Nystrom error", "Nystrom fit", "RK4 error", "RK4 fit", 'FontSize',11, "Location", "SOUTHEAST")

    myfig = gcf();
    exportgraphics(myfig , 'global_err_explicit_01.png', 'Resolution', 600);




   
end