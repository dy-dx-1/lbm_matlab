function [cd, cl] = aero_coeffs(f, obst_indices, dh, dt, calc_coeff)
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
    delta_P = 0; 
    for i=2:9
        summed_opps = f_reshaped(i,obst_indices) + f_reshaped(opp_cs(i),obst_indices);
        sum_obj_elements = sum(summed_opps); 
        delta_P = delta_P + sum_obj_elements*c(i,:); 
        disp(delta_P); 
    end    
    return 
    delta_P = delta_P*(dh^2); 
    F = delta_P/dt;  % in the shape of [Fx Fy] 
    % NOTE: CHECK AND CONFIRM UNIT MATCHES WITH CALC_COEFF!!!!!!!!!
    coeffs = (calc_coeff*F);  
    cd = coeffs(1); 
    cl = coeffs(2); 
end