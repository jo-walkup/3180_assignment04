%example for how to use compute_planetary_motion(...)
function usage_example_planet_motion()
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

    G_rate= @(t,v) gravity_rate_func(t, v,orbit_params);
    [t_list,X_list,h_avg, num_evals, failed_steps_rate] = explicit_RK_variable_step_integration(G_rate,t_range,V0,h_ref,DormandPrince,p,error_desired);

    
    subplot(2,1,1);

    plot(t_range, V_list(:, 1), 'g'); hold on
    plot(t_list,X_list(:,1),'k.');
 

    xlabel("Time (-)", FontSize=13);
    ylabel("X position(-)", FontSize=13);
    title("Comparison of Methods", FontSize=13);
    
    subplot(2,1,2); 

    plot(t_range, V_list(:, 2), 'g'); hold on
    plot(t_list,X_list(:,2),'k.');

    xlabel("Time (-)", FontSize=13);
    ylabel("Y position (-)", FontSize=13);

    mylegend = legend("Analytical", "Runge Kuntta", FontSize=11);
    set(mylegend, "location", "southoutside", "orientation", "horizontal")

    myfig = gcf();
    exportgraphics(myfig , 'planetary_motionxandy.png', 'Resolution', 600);
end