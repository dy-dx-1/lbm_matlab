function [f] = apply_meso_obs(f, u_in, b_obs_indices)
% Applies mesoscopic scale BCs for the general case of solid top and bottom
% walls, west constant speed inlet and east constant pressure outlet 
    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = pressure_bc(f,'east'); 
    %f = outlet_bc(f, 'east');  % produces no reverberations|!!!!
    f = inlet_bc(f, u_in, 'west'); % constant entry speed at the left 
    f = obstacle_bb(f, b_obs_indices); 
end