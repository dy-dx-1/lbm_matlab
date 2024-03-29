function show_velocity(show_vector_field, u, v, u_lb, sample_factor, x_sampled, y_sampled, rho, pressure_calc_coeff, a_cyl_indices) 
    % Displays only the velocity graph 
    uu = sqrt(u.^2+v.^2) / u_lb;
    uu(a_cyl_indices) = nan; 
    imagesc(flipud(uu));
    if show_vector_field==1
        % Arrow vector field 
        hold on;
        sample_u = flipud(u(1:sample_factor:end, 1:sample_factor:end));
        sample_v = flipud(v(1:sample_factor:end, 1:sample_factor:end));
        quiver(x_sampled, y_sampled, sample_u, sample_v, 'r'); 
        hold off; 
    end 
    title("Champ de vitesse normalisé par rapport à u_{lb}")
    xlabel("Noeuds en X") 
    ylabel("Noeuds en Y") 
    colorbar
    axis equal; 