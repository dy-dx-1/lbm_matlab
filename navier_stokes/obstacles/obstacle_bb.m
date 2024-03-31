function f = obstacle_bb(f, b_obs_indices, boundary_links)
    % D2Q9
    % Applies a wall bounceback BC to an obstacle's linear indices 
    % f: distribution function
    % b_obs_indices: linear indices of the boundary nodes of the obstacle. Since we're using halfway BB, 
    % the actual frontier is located half a step away, so these nodes are technically fluid nodes that are about to collide with the obstacle.
    
    opp_cs = [ 1,   4,  5,  2,  3,  8,   9,   6,   7]; % indices of opposite speeds for bb 
    f_reshaped = permute(f, [3 1 2]);
    [~, ny, nx] = size(f_reshaped); 
    for i=1:9
        f_reshaped(i, b_obs_indices) = f_reshaped(opp_cs(i), b_obs_indices);
    % Iterating through all boundary nodes defined by b_obst_indices
    %{
    for i=1:length(b_obs_indices)
        % converting linear indexation of obst_indices to 2D indexation to find corresponding links in boundary_links
        [row, col] = ind2sub([ny, nx], b_obs_indices(i));
        links = boundary_links{row, col}; % we now have a  1D array [] of all links connected to the boundary node
        for j=1:length(links)
            link = links(j); % This now represents a direction in the D2Q9 lattice where we want to apply BB
            f_reshaped(link, b_obs_indices(i)) = f_reshaped(opp_cs(link), b_obs_indices(i));
        end
    end
    %}
    f = permute(f_reshaped, [2 3 1]);
end