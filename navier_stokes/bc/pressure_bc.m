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
        
        c = zeros(9,2);
        c(1,:) = [0, 0];
        c(2,:) = [1, 0];
        c(3,:) = [0, 1];
        c(4,:) = [-1, 0];
        c(5,:) = [0, -1];
        c(6,:) = [1, 1];
        c(7,:) = [-1, 1];
        c(8,:) = [-1, -1];
        c(9,:) = [1, -1];
        
        % vitesse sur le mur, vector form 
        rho = sum(f,3);
        [rows, cols] = size(rho);

        u_global, v_global, rho_ = reconstruct_macro_all(f); 
        u_b = u_global(2:end-1, end-1); 
        v_b = v_global(2:end-1, end-1);
        u_vect_b = [u_b, v_b]; 

        u_b1 = u_global(2:end-1, end); 
        v_b1 = v_global(2:end-1, end);
        u_vect_b1 = [u_b1, v_b1];

        u_vect_w = u_vect_b + 0.5*(u_vect_b-u_vect_b1); 
        % (prendre tt donnees verticales de avant derniere col pr BB) 
        f(2:end-1, end-1, 4) = -f(2:end-1, end-1, 2) + (2*w_2*rho_boundary*(1+(((dot(c(2,:), u_vect_w)))/(2*(cs^4)))-((u_vect_w^2)/(2*(cs^2)))));
        f(2:end-1, end-1, 7) = -f(2:end-1, end-1, 6) + (2*w_6*rho_boundary*(1+(((dot(c(2,:), u_vect_w)))/(2*(cs^4)))-((u_vect_w^2)/(2*(cs^2)))));
        f(2:end-1, end-1, 8) = -f(2:end-1, end-1, 9) + (2*w_9*rho_boundary*(1+(((dot(c(2,:), u_vect_w)))/(2*(cs^4)))-((u_vect_w^2)/(2*(cs^2)))));

    end
    
    