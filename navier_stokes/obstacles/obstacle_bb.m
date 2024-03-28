function f = obstacle_bb(f, b_obs_indices)
    % D2Q9
    % Applies a wall bounceback BC to an obstacle's linear indices 
    % f: distribution function
    % b_obs_indices: linear indices of the boundary nodes of the obstacle. Since we're using halfway BB, 
    % the actual frontier is located half a step away, so these nodes are technically fluid nodes that are about to collide with the obstacle.
    
    opp_cs = [ 1,   4,  5,  2,  3,  8,   9,   6,   7]; % indices of opposite speeds for bb 
    f_reshaped = permute(f, [3 1 2]);
    for i=1:9
        f_reshaped(i, b_obs_indices) = f_reshaped(opp_cs(i), b_obs_indices);  
    end     
    f = permute(f_reshaped, [2 3 1]);
end