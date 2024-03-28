function links_matrix = find_boundary_links(ins_obs_matrix)
    % Iterates over the domain and finds the nodes defining the boundary. It then associates the connecting links 
    % to the boundary nodes. Needed for momentum exchange method. 
    % Note: the ins_obs_matrix MUST not include the actual boundaries 
    % (i.e. must not be the full_obs_matrix) because we use the nearest
    % 0s to find the connecting links to the solid nodes 1s 
    % This means that actual boundary must be 0 in this code for the proper
    % links to be associated. If the full matrix is used, the links will be
    % generated for nodes shifted 1 place to the exterior (that's where
    % the zeros will be) and aero_coeffs will just have an empty
    % links array when it tries to index it's boundary nodes. 
    [ny, nx] = size(ins_obs_matrix);
    links_matrix = cell(ny, nx);
    
    for i = 2:ny-1
        for j = 2:nx-1
            if ins_obs_matrix(i, j) == 0
                links = [];
                
                if ins_obs_matrix(i-1, j) == 1
                    links = [links, 3];
                end
                
                if ins_obs_matrix(i+1, j) == 1
                    links = [links, 5];
                end
                
                if ins_obs_matrix(i, j-1) == 1
                    links = [links, 4];
                end
                
                if ins_obs_matrix(i, j+1) == 1
                    links = [links, 2];
                end
                
                if ins_obs_matrix(i+1, j-1) == 1
                    links = [links, 8];
                end
                
                if ins_obs_matrix(i+1, j+1) == 1
                    links = [links, 9];
                end
                
                if ins_obs_matrix(i-1, j-1) == 1
                    links = [links, 7];
                end
                
                if ins_obs_matrix(i-1, j+1) == 1
                    links = [links, 6];
                end
                
                links_matrix{i, j} = links;
            end
        end
    end
end