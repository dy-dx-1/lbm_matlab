% A Lattice Boltzmann (single relaxation time) D2Q9 solver,
% with the Spalart Allmaras turbulence model, on a lid-driven cavity. 
% Cell centers (nodes) are placed on the boundaries. 
% Author: Robert Lee
% Email: rlee32@gatech.edu

% Modifications au code a  des fins de validation faites par Nicolas Sarabia-Benoit 
% Dans le cadre du cours MEC3900 - Projet Integrateur III, a  Polytechnique Montreal 

clear;close all;clc;

addpath basic
addpath bc
addpath turbulence
addpath verif_assets
addpath obstacles

%%% Base parameters.
Re = 100; % Nombre de Reynolds, a  commenter pour imposer viscosite cinematique 
tau = 0.809; 
rho0 = 1; 
total_time = 20; 
% Geometric parameters.
nx = 130; % total # of nodes used in the sim, 16641 corresponds to 129x129 for GHIA comparison 
duct_ratio = 2; % ratio for the rectangel such that Length = ratio*Height 
if mod(nx,duct_ratio) ~= 0
    error('nx must be divisible by duct_ratio to keep dx=dy');
end
cyl_size_ratio = 0.2; % Diam of cyl as a fraction of the duct height

%%% Derived parameters
ny = nx/duct_ratio; 
nu_lb = (tau-0.5)/3; % kinematic viscosity in lb units
dh = 1/(nx-1); % Spacing between nodes in the x and y direction 
dt = nu_lb*(dh^2)*Re;
u_lb = dt/dh; % lattice speed 
timesteps = round(total_time/dt); % total number of timesteps
omega = 1/tau; % relaxation parameter

%%% Simulation parameters.
% Setting up cylinder area 
[X,Y] = meshgrid(1:nx,1:ny);
cyl_rad_nodes = round(cyl_size_ratio*ny*0.5); % cyl radius expressed in nodes
x_cyl = round(nx/3); % X position of center of cyl 
y_cyl = round(ny/2); % Y position of center of cyl
cyl_matrix = generate_obstacle_matrix(X, Y, x_cyl, y_cyl, cyl_rad_nodes, 'circle');  % Matrix where 1 represents a cylinder node 
cyl_indices = find(cyl_matrix); % linear indexation of non zero elemetns, will be used to apply BB https://www.mathworks.com/help/matlab/ref/find.html

% Displaying info 
display_sim_info(nx, ny, timesteps, Re, dh, dt, nu_lb, tau);

% Initialize.
rho = rho0*ones(ny,nx);
u = zeros(ny,nx);
v = zeros(ny,nx);
u(2:end-1,1) = u_lb;

f = compute_feq(rho,u,v);
% Apply meso BCs.
f = wall_bc(f,'north');
f = wall_bc(f,'south');
f = pressure_bc(f,'east'); 
f = inlet_bc(f, u_lb, 'west'); % constant entry speed at the left 
f = obstacle_bb(f, cyl_indices); 

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
u(2:end-1,1) = u_lb;
v(2:end-1,1) = 0;
% East wall (outlet), constant pressure 
rho(2:end-1, end) = 1; % constant density 
v(2:end-1, end) = 0; % uniform flow out 


% Main loop.
disp(['Running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps
    if (mod(iter,round(timesteps/10))==0)
        disp(['Running ... ' num2str(iter/timesteps*100) '% completed']);
    end
    
    % Collision.
    f = collide_mrt(f, u, v, rho, omega);
    
    % Apply meso BCs.
    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = pressure_bc(f,'east'); 
    f = inlet_bc(f, u_lb, 'west'); % constant entry speed at the left 
    f = obstacle_bb(f, cyl_indices); 

    % Streaming.
    f = stream(f);
    
    % Apply meso BCs.
    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = pressure_bc(f,'east');  
    f = inlet_bc(f, u_lb, 'west'); % constant entry speed at the left 
    f = obstacle_bb(f, cyl_indices); 
    
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
    u(2:end-1,1) = u_lb;
    v(2:end-1,1) = 0;
    % East wall (outlet), constant pressure 
    rho(2:end-1, end) = 1; % constant density 
    v(2:end-1, end) = 0; % uniform flow out 
    % super instable u(2:end-1, end) = -1 + (1./rho(2:end-1, end)).*sum(f(2:end-1, end, [1,3,5]), 3) + 2*sum(f(2:end-1, end, [2,6,9]), 3);
    
    % Sanity check that the simulation is working, i.e checking that all non boundary nodes in the u matrix are not NaN
    if (any(isnan(u(2:end-1,2:end-1))))
        disp('!!!!!!!!!!!!!!!---- Error: NaNs in u matrix. Exiting ----!!!!!!!!!!!!!!!');
        return;
    end    
    % VISUALIZATION
    % Modified from Jonas Latt's cavity code on the Palabos website.
    if (mod(iter,10==0))
        uu = sqrt(u.^2+v.^2) / u_lb;
        uu(cyl_indices) = nan; 
        imagesc(flipud(uu));
%       imagesc(flipud(nut)); %%% turbulence viscosity 
%       imagesc(flipud(omega));
        % rectangle function is easiest to draw a circle, pos vector outlines lower left corner and height and width, curvature makes it a circle 
        rectangle('Position', [x_cyl-cyl_rad_nodes y_cyl-cyl_rad_nodes cyl_rad_nodes*2 cyl_rad_nodes*2], 'Curvature', [1 1], 'FaceColor', 'red')
        colorbar
        axis equal; 
        drawnow
    end
end
%{
    TODO: FORMAT APPROPIATELY SO THAT IT OUTPUTS 2 COLUMSN EVEN THOUGH NX!=NY
% Exctracting velocity data along the middle for validation with GHIA
u_center = flipud(extractRowOrColumn(u, 'col', round(nx/2)))/u_lb; % u along the vertical line at center 
v_center = flipud(extractRowOrColumn(v, 'row', round(ny/2)))/u_lb; % v along the horizontal line at center 
% displaying velocities in array form 
out_center_speeds = zeros(nodes, 2); 
out_center_speeds(:,1) = u_center; 
out_center_speeds(:,2) = v_center; 
disp(out_center_speeds); 
disp('**************Done!****************');
disp('Done!');
%} 