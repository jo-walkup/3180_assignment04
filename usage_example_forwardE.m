%example for how to use compute_planetary_motion(...)
function usage_example_forwardE()
     

    BT=struct();

    BT.A=[0];
    BT.B=[1];
    BT.C=[0];
    
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
    h_ref=0.1;

    G_rate= @(t,v) gravity_rate_func(t, v,orbit_params);
    
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(G_rate,t_range,V0,h_ref,BT);

    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(V_list(:,1),V_list(:,2),'k');

    plot(X_list(:,1),X_list(:,2),'b');

    legend("Sun")
end