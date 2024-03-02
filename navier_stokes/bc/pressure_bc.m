function f = pressure_bc(f, side)
    % D2Q9
    % Applies a pressure BC to the distirbution functions in 2d matrix form.
    % Assumes boundary nodes are placed half a cell away. NEBB wet node scheme not usable because it is not appropiate for high Re. 
    % NOTE: i from 1 to 9 (and not 0 to 8) 
    
    if strcmp(side, 'east') % East outlet.   
        rho_boundary = 1; % densite constante, ecoulement incompressible 
        cs = sqrt(1/3); 

        % poids et directions
        w_2 = 1/9;
        w_6 = 1/36; 
        w_9 = 1/36;  
        
        c_2 = [1; 0];
        c_6 = [1; 1];
        c_9 = [1; -1];
        
        % vitesse sur le mur, vector form 
        [u, v, rho] = reconstruct_macro_all(f); 
        rho_boundary = rho(2:end-1, end); 
        disp(rho_boundary); 
        u_b = u(2:end-1, end); %127x1 
        v_b = v(2:end-1, end); %127x1 
        u_vect_b = [u_b, v_b]; 
        
        u_b1 = u(2:end-1, end-1); % ub and ub1 refer to the boundary node and the next interior node following the inward normal vector of the boundary.
        v_b1 = v(2:end-1, end-1);
        u_vect_b1 = [u_b1, v_b1];

        u_vect_w = u_vect_b + 0.5*(u_vect_b-u_vect_b1); % vect [nodes-2 x 2] of approx speeds u, v along boundary 

        % (prendre tt donnees verticales de avant derniere col pr BB) 
        f(2:end-1, end, 4) = -f(2:end-1, end, 2) + (2*w_2*rho_boundary.*(1+(((u_vect_w*c_2).^2)/(2*(cs^4)))-((vecnorm(u_vect_w, 2, 2).^2)/(2*(cs^2)))));
        f(2:end-1, end, 7) = -f(2:end-1, end, 9) + (2*w_9*rho_boundary.*(1+(((u_vect_w*c_9).^2)/(2*(cs^4)))-((vecnorm(u_vect_w, 2, 2).^2)/(2*(cs^2)))));
        f(2:end-1, end, 8) = -f(2:end-1, end, 6) + (2*w_6*rho_boundary.*(1+(((u_vect_w*c_6).^2)/(2*(cs^4)))-((vecnorm(u_vect_w, 2, 2).^2)/(2*(cs^2)))));

    end