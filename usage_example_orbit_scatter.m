%example for how to use compute_planetary_motion(...)
function usage_example_orbit_scatter()
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    clf();
    clear all;

    BT_fe=struct();

    BT_fe.A=[0];
    BT_fe.B=[1];
    BT_fe.C=[0];

    BT_rk4=struct();

    BT_rk4.A=[0, 0, 0, 0; 1/2, 0, 0, 0; 0, 1/2, 0, 0; 0, 0, 1, 0];
    BT_rk4.B=[1/6; 1/3; 1/3; 1/6];
    BT_rk4.C=[0; 0.5; 0.5; 1];

    BT_em=struct();

    BT_em.A=[0, 0; 0.5, 0];
    BT_em.B=[0;1];
    BT_em.C=[0; 0.5];
    
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
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

    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);

    h_ref=0.0005;
    p=5.5;
    error_desired=0.0001;

    hs = logspace(-3,-1,100);

    G_rate= @(t,v) gravity_rate_func(t, v,orbit_params);

    h_lists = zeros(size(hs));
    distances = zeros(size(hs));


    % for n = 1:length(hs)
        % h = hs(n);
        [t_list,X_list,h_avg, num_evals, failed_steps_rate] = explicit_RK_variable_step_integration(G_rate,t_range,V0,h_ref,DormandPrince,p,error_desired);
        h_list = diff(t_list);

        x_vals=X_list(:,1);
        y_vals=X_list(:,2);
        distance=sqrt(x_vals.^2+y_vals.^2);

        % h_lists(n) = h_list
        % distances(n) = distance;

    % end
    distance
    h_list

    loglog(distance(1:59),h_list, "b.", MarkerSize=20);

    myfig = gcf();
    exportgraphics(myfig , 'planetary_motion_scatter.png', 'Resolution', 600);
end