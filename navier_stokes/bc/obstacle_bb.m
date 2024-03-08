function f = obstacle_bb(f, obst_indices)
    % D2Q9
    % Applies a wall bounceback BC to an obstacle's linear indices 
    opp_cs = [ 1,   4,  5,  2,  3,  8,   9,   6,   7]; % indices of opposite speeds for bb 
    for i=1:9
        f(obst_indices, i) = f(obst_indices, opp_cs(i));
    end     
    end