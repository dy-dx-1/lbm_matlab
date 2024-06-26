function display_sim_info(nx, ny, timesteps, Re, dh, dt, nu_lb, tau, cyl_rad_nodes, cyl_size_ratio)
    % Displays relevant simulation information 
    disp("-------- Simulation started --------"); 
    disp("*** Simulation parameters ***"); 
    disp(strcat("Reynolds number: ", num2str(Re))); 
    disp(strcat("Relaxation time tau: ", num2str(tau)));
    disp(strcat("Nodes in the X direction: ", num2str(nx)));
    disp(strcat("Nodes in the Y direction: ", num2str(ny)));
    disp(strcat("Total time simulated: ", num2str(round(dt*timesteps))));
    disp(strcat("Dim radius of cylinder dx*(rad_cyl_nodes-1): ", num2str((cyl_rad_nodes-1)*dh)));
    disp(strcat("Cyl size ratio: D_cyl = ", num2str(cyl_size_ratio), "*H || H = ", num2str(1/cyl_size_ratio), "*D_cyl")); 
    disp("*** Derived parameters ***"); 
    disp(strcat("Spatial step dx=dy: ", num2str(dh))); 
    disp(strcat("Timestep dt: ", num2str(dt))); 
    disp(strcat("LB viscosity: ", num2str(nu_lb))); 
     
    ratio_relax_dt = tau/dt; 
    if (ratio_relax_dt < 1)
        disp(strcat("WARNING! tau/dt ratio lower than 1, ref stab checks for BGK unavailable. tau/dt = ", num2str(ratio_relax_dt)))
    end 
    if (tau > 5 && tau < 0.5) 
        disp("WARNING!! Relaxation time tau not stable"); 
    end 
end