function display_sim_info(L, U, nodes, timesteps, Re, dx, dt, u_lb, tau)
    % Displays relevant simulation information 
    disp("--------Simulation started--------"); 
    disp("***Simulation parameters***"); 
    disp(["Cavity length: " num2str(L)]); 
    disp(["Lid velocity: " num2str(U)]); 
    disp(["Reynolds number: " num2str(Re)]); 
    disp(["Number of nodes: " num2str(nodes)]); 
    disp(["Number of timesteps: " num2str(timesteps)]); 
    disp(["***Derived parameters***"]); 
    disp(["Spatial step dx: " num2str(dx)]); 
    disp(["Timestep dt: " num2str(dt)]); 
    disp(["Total time simulated: " num2str(dt*timesteps)]);
    disp(["----"]); 
    disp(["Lattice speed: " num2str(u_lb)]); 
    disp(["Relaxation time: " num2str(tau)]); 
    disp(["***Stability checks***"]); 
    ratio_relax_dt = tau/dt; 
    if (ratio_relax_dt<1)
        disp(strcat("WARNING! tau/dt ratio lower than 1, ref stab checks for BGK unavailable. tau/dt = ", num2str(ratio_relax_dt)))
    end 
    cond_stab_bgk = sqrt(2/3)*(dx/dt); % juste pour avoir un idee par rapport a  BGK
    disp(strcat("BGK stability condition (should be >U_p): ", num2str(cond_stab_bgk))); 
    if (U>cond_stab_bgk)
        disp("WARNING!!! BGK stab condition not respected"); 
    end 
    if (tau>5 && tau<0.5) 
        disp("WARNING!! Relaxation time tau not stable"); 
    end 