function show_pressure(show_vector_field, u, v, u_lb, sample_factor, x_sampled, y_sampled, pressure_calc_coeff, a_cyl_indices) 
    % Getting pressure field & displaying it 
    p = rho.*pressure_calc_coeff;      
    contourf(p, 50);
    title("Champ de pression")
    xlabel("Noeuds en X") 
    ylabel("Noeuds en Y") 
    colorbar
    axis equal; 