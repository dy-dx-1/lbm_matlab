function obj_matrix = generate_obstacle_matrix(X_mesh, Y_mesh, center_x, center_y, caracteristic_length, type) 
    if type == "circle"
    % Caracteristic length is expected to be the radius of the circle
    obj_matrix = ((X_mesh-center_x).^2 + (Y_mesh-center_y).^2) <= (caracteristic_length.^2);
    end 
    if type == "rectangle"
      % caracteristic length is expected to be the length of a side of rectangle
      obj_matrix = (X_mesh >= center_x - caracteristic_length) & (X_mesh <= center_x + caracteristic_length) & (Y_mesh >= center_y - caracteristic_length) & (Y_mesh <= center_y + caracteristic_length);
    end 