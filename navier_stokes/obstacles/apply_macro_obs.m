function [u, v, rho] = apply_macro_obs(f, u_in)
    addpath basic
    % Applies macroscopic scale BCs for the general case of solid top and bottom
    % walls, west constant speed inlet and east constant pressure outlet
    
    % Determine macro variables
    [u,v,rho] = reconstruct_macro_all(f);
    % Apply Macro BCs 
    % North wall 
    u(end,:) = 0; 
    v(end,:) = 0;
    % South wall 
    u(1,:) = 0;
    v(1,:) = 0;
    % West wall (inlet)
    u(2:end-1,1) = u_in;
    v(2:end-1,1) = 0;
    % East wall (outlet), constant pressure 
    rho(2:end-1, end) = 1; % constant density 
    v(2:end-1, end) = 0; % uniform flow out 
    % super instable u(2:end-1, end) = -1 + (1./rho(2:end-1, end)).*sum(f(2:end-1, end, [1,3,5]), 3) + 2*sum(f(2:end-1, end, [2,6,9]), 3);
 
end


