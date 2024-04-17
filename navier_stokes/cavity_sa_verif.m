% UNDER CONSTRUCTION

% A Lattice Boltzmann Multiple Relaxation Time D2Q9 solver,
% with a Spalart Allmaras turbulence model, on a lid-driven cavity.
% This features a non-lattice-aligned wall! 
% Cell centers (nodes) are placed on the boundaries. 
% Author: Robert Lee
% Email: rlee32@gatech.edu

clear;close all;clc;

addpath basic
addpath bc
addpath verif_assets

% Algorithm steps:
% Initialize meso (f)
% Apply meso BCs
% Determine macro variables and apply macro BCs
% Loop:
%   Collide
%   Apply meso BCs
%   Stream
%   Apply meso BCs?
%   Determine macro variables and apply macro BCs

% Physical parameters.
Re = 7500; 
rho0 = 1;
% Discrete/numerical parameters.
nodes = 229;
total_time = 200; % in seconds 

tau = 0.535; % relaxation time 
nu_lb = (tau-0.5)/3; % kinematic viscosity in lb units

% Derived discrete parameters.
dh = 1/(nodes-1); 
dt = nu_lb*(dh^2)*Re;

omega = 1 / tau;
u_lb = dt / dh;
timesteps = round(total_time/dt); 


% Initialize.
f = ones(nodes,nodes,9);
nutilde = 1.2e-3*ones(nodes,nodes);
% Apply meso BCs.
f = moving_wall_bc(f,'north',u_lb);
f = wall_bc(f,'south');
f = wall_bc(f,'east');
f = wall_bc(f,'west');
% Determine macro variables and apply macro BCs
[u,v,rho] = reconstruct_macro_all(f);
u(end,2:end-1) = u_lb;
v(end,2:end-1) = 0;
u(1,:) = 0;
v(1,:) = 0;
u(:,1) = 0;
v(:,1) = 0;
u(:,end) = 0;
v(:,end) = 0;
nutilde(1,:) = 0;
nutilde(end,:) = 0;
nutilde(:,1) = 0;
nutilde(:,end) = 0;

% Main loop.
disp(['Running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps
    if (mod(iter,round(timesteps/10))==0)
        disp(['Running ... ' num2str(iter/timesteps*100) '% completed']);
    end
    
    % Collision.
    f = collide_mrt(f, u, v, rho, omega);
    
    % Apply meso BCs.
    f = moving_wall_bc(f,'north',u_lb);
    f = wall_bc(f,'south');
    f = wall_bc(f,'east');
    f = wall_bc(f,'west');

    % Streaming.
    f = stream(f);
    
    % Apply meso BCs.
    f = moving_wall_bc(f,'north',u_lb);
    f = wall_bc(f,'south');
    f = wall_bc(f,'east');
    f = wall_bc(f,'west');
    
    % Determine macro variables and apply macro BCs
    [u,v,rho] = reconstruct_macro_all(f);
    u(end,2:end-1) = u_lb;
    v(end,2:end-1) = 0;
    u(1,:) = 0;
    v(1,:) = 0;
    u(:,1) = 0;
    v(:,1) = 0;
    u(:,end) = 0;
    v(:,end) = 0;

    average_density = mean(rho, 'all');
    fileID = fopen('C:/Users/Nicolas/Downloads/average_density.txt', 'a');
    fprintf(fileID, '%f\n', average_density); 
    fclose(fileID); 
    disp(average_density); % displaying average density to check for conservation
    % Sanity check that the simulation is working, i.e checking that all non boundary nodes in the u matrix are not NaN
    if (any(isnan(u(2:end-1,2:end-1))))
        disp('!!!!!!!!!!!!!!!---- Error: NaNs in u matrix. Exiting ----!!!!!!!!!!!!!!!');
        return;
    end    
    % VISUALIZATION
    % Modified from Jonas Latt's cavity code on the Palabos website.
    %if (mod(iter,10)==0) uncomment to only visualise every other iter
        subplot(2,1,1); 
        uu = sqrt(u.^2+v.^2) / u_lb; 
        imagesc(flipud(uu));
        colorbar
        axis equal; 
        
        subplot(2,1,2); 
        % Getting pressure field & displaying it 
        p = rho.*(1/3)*((dh/dt))^2;      
        contourf(p, 50);
        colorbar
        axis equal; 
        drawnow
    %end
end
% Exctracting velocity data along the middle for validation with GHIA
u_center = flipud(extractRowOrColumn(u, 'col', round(nodes/2)))/u_lb; % u along the vertical line at center 
v_center = flipud(extractRowOrColumn(v, 'row', round(nodes/2)))/u_lb; % v along the horizontal line at center 
% displaying velocities in array form 
out_center_speeds = zeros(nodes, 2); 
out_center_speeds(:,1) = u_center; 
out_center_speeds(:,2) = v_center; 
disp(out_center_speeds); 
disp('**************Done!****************');
disp('Done!');