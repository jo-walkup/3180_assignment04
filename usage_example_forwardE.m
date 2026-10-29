%example for how to use compute_planetary_motion(...)
function usage_example_forwardE()
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    clf();
    clear all;

    % Ralston
    BT_fe=struct();

    BT_fe.A=[0, 0; 2/3, 0];
    BT_fe.B=[1/4; 3/4];
    BT_fe.C=[0; 2/3];

    BT_rk4=struct();

    BT_rk4.A=[0, 0, 0, 0; 1/2, 0, 0, 0; 0, 1/2, 0, 0; 0, 0, 1, 0];
    BT_rk4.B=[1/6; 1/3; 1/3; 1/6];
    BT_rk4.C=[0; 0.5; 0.5; 1];

    % Heun
    BT_em=struct();

    BT_em.A=[0, 0; 1, 0];
    BT_em.B=[1/2;1/2];
    BT_em.C=[0; 1];
    
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
    h_ref=0.0005;

    G_rate= @(t,v) gravity_rate_func(t, v,orbit_params);
    
    [t_list_fe,X_list_fe,h_avg, num_evals] = explicit_RK_fixed_step_integration(G_rate,t_range,V0,h_ref,BT_fe);
    [t_list_rk4,X_list_rk4,h_avg, num_evals] = explicit_RK_fixed_step_integration(G_rate,t_range,V0,h_ref,BT_rk4);
    [t_list_em,X_list_em,h_avg, num_evals] = explicit_RK_fixed_step_integration(G_rate,t_range,V0,h_ref,BT_em);

    % axis equal; axis square;
    % axis([-20,20,-20,20])
    % hold on
    % plot(0,0,'ro','markerfacecolor','r','markersize',5);
    % plot(V_list(:,1),V_list(:,2),'k'); hold on
    % 
    % plot(X_list_fe(:,1),X_list_fe(:,2),'b', 'markersize',20);
    % plot(X_list_rk4(:,1),X_list_rk4(:,2),'g', 'markersize',20);
    % plot(X_list_em(:,1),X_list_em(:,2),'m', 'markersize',20);
    % 
    % legend("Sun", "Analytical", "Forward Euler", "RK4", "Explicit Midpoint", FontSize=11)
    % title("Comparison of Methods", FontSize=13);
    % xlabel("X position of planet (-)", FontSize=13);
    % ylabel("Y position of planet (-)", FontSize=13);

    subplot(2,1,1);
    plot(t_range, V_list(:, 1), 'g'); hold on
    plot(t_list_fe,X_list_fe(:,1),'k');
    plot(t_list_rk4,X_list_rk4(:,1),'b');
    plot(t_list_em,X_list_em(:,1),'m--');
    xlabel("Time (-)", FontSize=13);
    ylabel("X position(-)", FontSize=13);
    title("Integration of a Planetary System via Several Runge-Kutta Methods", FontSize=13);
    
    subplot(2,1,2); 
    plot(t_range, V_list(:, 2), 'g'); hold on
    plot(t_list_fe,X_list_fe(:,2),'k');
    plot(t_list_rk4,X_list_rk4(:,2),'b');
    plot(t_list_em,X_list_em(:,2),'m--');
    xlabel("Time (-)", FontSize=13);
    ylabel("Y position (-)", FontSize=13);

    mylegend = legend("Analytical", "Ralston", "RK4", "Heun", FontSize=11);
    set(mylegend, "location", "southoutside", "orientation", "horizontal")

    myfig = gcf();
    exportgraphics(myfig , 'planet_path_all_methods_subplot.png', 'Resolution', 600);
end