function [u, v, rho, f] = run_LBM_loop(f, u, v, rho, omega, u_lb, b_cyl_indices, i_cyl_indices, a_cyl_indices, boundary_links, dh, dt, pressure_calc_coeff, p_inf, p_divider, timesteps, calc_coeff, sample_factor, x_cyl, y_cyl, cyl_rad_nodes, x_sampled, y_sampled, update_every_iter, vis, show_vector_field)
    % Main loop for LBM simulation
    % Multiple loops were prefered to avoid multiple if statements inside the loop (since we're running so many iterations).
    aero_coeffs = vis{1};
    density = vis{2};
    display_graphs = vis{3};
    cylinder = vis{4};
    
    if update_every_iter==1
        for iter = 1:timesteps    
            % Collision
            f = collide_mrt(f, u, v, rho, omega);
            
            % Apply meso BCs
            f = apply_meso_obs(f, u_lb, a_cyl_indices, boundary_links); 
        
            % Calculation of drag and lift coefficient 
            if (mod(iter, 10)) == 0  
                aero_coeffs(f, b_cyl_indices, boundary_links, dh, dt, calc_coeff, rho, pressure_calc_coeff, p_inf, p_divider); % calling the function to calculate the aero coeffs
            end
            % Streaming
            f = stream(f);    
            
            % Apply meso BCs
            f = apply_meso_obs(f, u_lb, a_cyl_indices, boundary_links); 
            
            % Apply macro variables
            [u,v,rho] = apply_macro_obs(f, u_lb, i_cyl_indices); 
            
            % Density check 
            density(rho); % calling the function to display avg density 
        
            % Sanity check that the simulation is working, i.e checking that all non boundary nodes in the u matrix are not NaN
            if (any(isnan(u(2:end-1,2:end-1))))
                disp('!!!!!!!!!!!!!!!---- Error: NaNs in u matrix. Exiting ----!!!!!!!!!!!!!!!');
                return;
            end    
            
            % VISUALIZATION & progress tracking 
            if (mod(iter,round(timesteps/10))==0)
                disp(['Running ... ' round(num2str(iter/timesteps*100)) '% completed']);
            end
            
            % calling the function to display the velocity||pressure||both field, putting all params even though they may not be used depending on the function 
            display_graphs(show_vector_field, u, v, u_lb, sample_factor, x_sampled, y_sampled, rho, pressure_calc_coeff, a_cyl_indices) 
            cylinder(x_cyl, y_cyl, cyl_rad_nodes); % calling the function to display the shape of the obstacle        
        
            drawnow
        end
    end 

    if update_every_iter==0
        for iter = 1:timesteps    
            % Collision
            f = collide_mrt(f, u, v, rho, omega);
            
            % Apply meso BCs
            f = apply_meso_obs(f, u_lb, a_cyl_indices, boundary_links); 
        
            % Calculation of drag and lift coefficient 
            if (mod(iter, 10)) == 0  
                aero_coeffs(f, b_cyl_indices, boundary_links, dh, dt, calc_coeff, rho, pressure_calc_coeff, p_inf, p_divider); % calling the function to calculate the aero coeffs
            end 
            % Streaming
            f = stream(f);    
            
            % Apply meso BCs
            f = apply_meso_obs(f, u_lb, a_cyl_indices, boundary_links); 
            
            % Apply macro variables
            [u,v,rho] = apply_macro_obs(f, u_lb, i_cyl_indices); 
            
            % Density check 
            density(rho); % calling the function to display avg density   
            
            % VISUALIZATION & progress tracking 
            if (mod(iter,round(timesteps/10))==0)
                disp(['Running ... ' round(num2str(iter/timesteps*100)) '% completed']);
                % Sanity check that the simulation is working, i.e checking that all non boundary nodes in the u matrix are not NaN
                if (any(isnan(u(2:end-1,2:end-1))))
                    disp('!!!!!!!!!!!!!!!---- Error: NaNs in u matrix. Exiting ----!!!!!!!!!!!!!!!');
                    return;
                end  
                
                % calling the function to display the velocity||pressure||both field, putting all params even though they may not be used depending on the function 
                display_graphs(show_vector_field, u, v, u_lb, sample_factor, x_sampled, y_sampled, rho, pressure_calc_coeff, a_cyl_indices) 
                cylinder(x_cyl, y_cyl, cyl_rad_nodes); % calling the function to display the shape of the obstacle        
            
                drawnow
            end 
        end
    end 