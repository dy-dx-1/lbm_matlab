function f = obstacle_bb(f, b_obs_indices, boundary_links)
    % D2Q9
    % Applies a wall bounceback BC to an obstacle's linear indices 
    % f: distribution function
    % b_obs_indices: linear indices of the boundary nodes of the obstacle. Since we're using halfway BB, 
    % the actual frontier is located half a step away, so these nodes are technically fluid nodes that are about to collide with the obstacle.
    
    opp_cs = [ 1,   4,  5,  2,  3,  8,   9,   6,   7]; % indices of opposite speeds for bb 
    [ny, nx, ~] = size(f); 
    for j=1:length(b_obs_indices)
        [row, col] = ind2sub([ny, nx], b_obs_indices(j));
        for i=1:9
            f(row, col, i) = f(row, col, opp_cs(i));
        end
    end 
end