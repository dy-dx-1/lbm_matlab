function display_sim_info(nx, ny, timesteps, Re, dh, dt, nu_lb, tau)
    % Displays relevant simulation information 
    disp("--------Simulation started--------"); 
    disp("***Simulation parameters***"); 
    disp(["Reynolds number: " num2str(Re)]); 
    disp(["Relaxation time tau: " num2str(tau)]);
    disp(["Nodes in the X direction: " num2str(nx)]);
    disp(["Nodes in the Y direction: " num2str(ny)]);
    disp(["Total time simulated: " num2str(round(dt*timesteps))]);
    disp(["***Derived parameters***"]); 
    disp(["Spatial step dx=dy: " num2str(dh)]); 
    disp(["Timestep dt: " num2str(dt)]); 
    disp(["LB viscosity: " num2str(nu_lb)]); 
     
    ratio_relax_dt = tau/dt; 
    if (ratio_relax_dt<1)
        disp(strcat("WARNING! tau/dt ratio lower than 1, ref stab checks for BGK unavailable. tau/dt = ", num2str(ratio_relax_dt)))
    end 
    if (tau>5 && tau<0.5) 
        disp("WARNING!! Relaxation time tau not stable"); 
    end 