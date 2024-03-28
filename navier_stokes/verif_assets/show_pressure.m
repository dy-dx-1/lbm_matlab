function show_pressure(show_vector_field, u, v, u_lb, sample_factor, x_sampled, y_sampled, rho, dh, dt, a_cyl_indices) 
    % Getting pressure field & displaying it 
    p = rho.*(1/3)*((dh/dt))^2;      
    contourf(p, 50);
    title("Champ de pression")
    xlabel("Noeuds en X") 
    ylabel("Noeuds en Y") 
    colorbar
    axis equal; 