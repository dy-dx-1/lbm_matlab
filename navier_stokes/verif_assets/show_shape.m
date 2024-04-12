function show_shape(x_cyl, y_cyl, cyl_rad_nodes)
    % rectangle function is easiest to draw a circle, pos vector outlines lower left corner and height and width, curvature makes it a circle 
    rectangle('Position', [x_cyl-cyl_rad_nodes y_cyl-cyl_rad_nodes-1 cyl_rad_nodes*2 cyl_rad_nodes*2], 'Curvature', [1 1], 'FaceColor', 'white')
