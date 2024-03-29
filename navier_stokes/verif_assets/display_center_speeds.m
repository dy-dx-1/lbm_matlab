function display_center_speeds(u, v, u_lb, nx, ny)
    % Extracts velocity data along the middle for validation with GHIA
    u_center = flipud(extractRowOrColumn(u, 'col', round(nx/2)))/u_lb; % u along the vertical line at center 
    v_center = flipud(extractRowOrColumn(v, 'row', round(ny/2)))/u_lb; % v along the horizontal line at center 
    % displaying velocities in array form 
    disp("FINAL SPEEDS ------------------");
    if nx==ny
        % Display side to side in a nodesx2 matrix
        out_center_speeds = zeros(nx, 2); 
        out_center_speeds(:,1) = u_center; 
        out_center_speeds(:,2) = v_center; 
        disp(out_center_speeds); 
    else 
        % Then display separately 
        disp("u along vertical line:");
        disp(u_center);
        disp("v along horizontal line:");
        disp(v_center);      
    end 
end 