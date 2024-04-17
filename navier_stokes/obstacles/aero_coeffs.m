function [cd, cl, cp] = aero_coeffs(f, obst_indices, boundary_links, dh, dt, calc_coeff, rho, pressure_calc_coeff, p_inf, p_divider)
    % D2Q9
    % Calculates the drag coefficient of the object defined by obst_indices
    % calc_coeff is a constant such that cd, cl = vect(F)*calc_coeff
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
    opp_cs = [1,   4,  5,  2,  3,  8,   9,   6,   7]; % indices of opposite speeds for bb 
    f_reshaped = permute(f, [3 1 2]);
    [~, ny, nx] = size(f_reshaped); 
    delta_P = [0, 0]; 
    % Iterating through all boundary nodes defined by obst_indices
    for i=1:length(obst_indices)
        % converting linear indexation of obst_indices to 2D indexation to find corresponding links in boundary_links
        [row, col] = ind2sub([ny, nx], obst_indices(i));
        links = boundary_links{row, col}; % we now have a  1D array [] of all links connected to the boundary node
        for j=1:length(links)
            link = links(j); % This now represents a direction in the D2Q9 lattice where we want to calculate momentum exchange
            % calculating momentum exhange for a singular link and projecting on the velocity vector
            summed_opps = (f_reshaped(link, obst_indices(i)) + f_reshaped(opp_cs(link), obst_indices(i)));
            delta_P = delta_P + summed_opps*c(link,:);
        end
    end
    F = delta_P; % [Fx, Fy] 
    %{
    Dimensionnal approach would be: 
    %delta_P = delta_P.*(dh^2); 
    %F = delta_P./dt;  % in the shape of [Fx Fy] 
    Here we use adim dx=dt=1 because Cd calc_coeff is defined with adim
    units in the main program 
    %}
    coeffs = (calc_coeff*F);  
    cd = coeffs(1); 
    cl = coeffs(2);
    %% Pressure coefficient 
    pressure_domain = rho.*pressure_calc_coeff; 
    pressure_boundary = pressure_domain(obst_indices); 
    cp = (pressure_boundary-p_inf)/(p_divider); 
end