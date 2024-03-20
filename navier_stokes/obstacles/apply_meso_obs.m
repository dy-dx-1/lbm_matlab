function [f] = apply_meso_obs(f, u_in, obstacle_indices)
    addpath bc
% Applies mesoscopic scale BCs for the general case of solid top and bottom
% walls, west constant speed inlet and east constant pressure outlet 
    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = pressure_bc(f,'east');  
    f = inlet_bc(f, u_in, 'west'); % constant entry speed at the left 
    f = obstacle_bb(f, obstacle_indices); 
end

