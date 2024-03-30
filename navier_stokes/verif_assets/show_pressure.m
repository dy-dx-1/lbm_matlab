function show_pressure(show_vector_field, u, v, u_lb, sample_factor, x_sampled, y_sampled, rho, pressure_calc_coeff, a_cyl_indices) 
    % Getting pressure field & displaying it 
    p = rho.*pressure_calc_coeff;      
    contourf(p, 50);
    % colorbar 