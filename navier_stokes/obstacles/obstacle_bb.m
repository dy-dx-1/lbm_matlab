function f = obstacle_bb(f, obst_indices)
    % D2Q9
    % Applies a wall bounceback BC to an obstacle's linear indices 
    opp_cs = [ 1,   4,  5,  2,  3,  8,   9,   6,   7]; % indices of opposite speeds for bb 
    f_reshaped = permute(f, [3 1 2]);
    for i=1:9
        f_reshaped(i, obst_indices) = f_reshaped(opp_cs(i), obst_indices);
    end     
    f = permute(f_reshaped, [2 3 1]);
end