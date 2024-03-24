function obj_matrix = mark_boundary_nodes(obj_matrix)
    % Initialize a matrix to store the boundary nodes
    boundary_nodes = zeros(size(obj_matrix));

    % Define the dimensions of the matrix
    [rows, cols] = size(obj_matrix);

    % Iterate over each point in the matrix
    for i = 2:rows-1
        for j = 2:cols-1
            % If the current point is inside the obstacle
            if obj_matrix(i, j) == 1
                % Check if any neighboring points are outside the obstacle
                if obj_matrix(i-1, j) == 0 || obj_matrix(i+1, j) == 0 || ...
                   obj_matrix(i, j-1) == 0 || obj_matrix(i, j+1) == 0
                    % Mark the current point as a boundary node
                    boundary_nodes(i, j) = 1;
                end
            end
        end
    end

    % Return the matrix containing only the boundary nodes
    obj_matrix = boundary_nodes;
end