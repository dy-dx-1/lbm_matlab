% Code simulating turbulent flow around a cylinder using the Lattice Boltzmann Method. 
% This code is based on the code in "cavity_sa.m" by Robert Lee (rlee32@gatech.edu)
% and modified by Nicolas Sarabia-Benoit (nicolas.sarabia-benoit@polymtl.ca) in the context
% of the course MEC3900 at Polytechnique Montreal.

clear;close all;clc;

addpath basic
addpath bc
addpath turbulence
addpath verif_assets
addpath obstacles

%%% Base parameters
Re = 40;          % Reynolds number
tau = 0.809;      % Relaxation time
total_time = 20;  % Total simulation time
rho0 = 1;         % Initial density (adim)

%%% Geometric parameters.
nx = 135;             % # of nodes in the x direction
duct_ratio = 3;       % ratio for the rectangular domain such that Length = ratio*Height 
cyl_size_ratio = 0.2; % Diam of cyl as a fraction of the duct height
if mod(nx,duct_ratio) ~= 0
    error('nx must be divisible by duct_ratio to keep dx=dy');
end

%%% Derived simulation parameters
ny = nx/duct_ratio; 
nu_lb = (tau-0.5)/3;  % kinematic viscosity in lb units
dh = 1/(nx-1);        % spacing between nodes in the x and y direction 
dt = nu_lb*(dh^2)*Re; % time step
u_lb = dt/dh;         % lattice speed 
timesteps = round(total_time/dt); % total number of timesteps
omega = 1/tau;                    % relaxation parameter

%%% Derived cylinder parameters
[X,Y] = meshgrid(1:nx,1:ny);
x_cyl = round(nx/3);                          % X position of center of cyl 
y_cyl = round(ny/2);                          % Y position of center of cyl
cyl_rad_nodes = round(cyl_size_ratio*ny*0.5); % cyl radius expressed in nodes
cyl_matrix = generate_obstacle_matrix(X, Y, x_cyl, y_cyl, cyl_rad_nodes, 'circle');  % Matrix where 1 represents a cylinder node 
boundary_cyl_matrix = mark_boundary_nodes(cyl_matrix);                               % Matrix where 1 represents ONLY the boundary of the cylinder
boundary_links = find_boundary_links(cyl_matrix-boundary_cyl_matrix);                % Getting the velocity links for each boundary node, used for momentum exchange
% Prepping linear indexes of non zero elements, will be used to easily index on cylinder https://www.mathworks.com/help/matlab/ref/find.htm
a_cyl_indices = find(cyl_matrix); % boundary and inside 
b_cyl_indices = find(boundary_cyl_matrix);  % only boundary
i_cyl_indices = find(cyl_matrix-boundary_cyl_matrix);  % only inside without boundary

%%% Displaying simulation info 
display_sim_info(nx, ny, timesteps, Re, dh, dt, nu_lb, tau, cyl_rad_nodes, cyl_size_ratio);

%%% Calculations for lift and drag coeff (see "obstacles/aero_coeffs.m") 
calc_coeff = (1/(rho0*(u_lb^2)*cyl_rad_nodes)); % cd, cl = vect(F)*calc_coeff

%%% Prepping for visualization (to avoid loop calculations)
% Velocity vector field (quiver)
sample_factor = 3;             % To have a less dense quiver field
[x_sampled, y_sampled] = meshgrid(1:sample_factor:nx, 1:sample_factor:ny); 
x_sampled = flipud(x_sampled); 
y_sampled = flipud(y_sampled); % flipped to match the image display
% Streamlines
N = 75; % number of seed locations ; used to display streamlines 
xstart = nx*rand(N,1); 
ystart = ny*rand(N,1);

%%% Visualization choices 
% This is a cell array that contains the visualization options for the simulation
% It needs to respect the following format: {aero_coeffs, density, velocity||pressure||velocity_pressure, show_shape}
% To not call a function use @placeholder
% Functions are defined in /verif_assets/
vis = {@placeholder, @placeholder, @show_velocity, @placeholder};
update_every_iter = 1; % 1 to update every iteration, 0 to update every 10% of the simulation, use 0 for large simulations
show_vector_field = 0; % 1 to show velocity quiver, 0 to not 

%%% Simulation ---------------------------------------------------------------------------------------------------------------------
% Initialization
rho = rho0*ones(ny,nx);
u = zeros(ny,nx);
v = zeros(ny,nx);
u(2:end-1,1) = u_lb;

f = compute_feq(rho,u,v);
f = apply_meso_obs(f, u_lb, b_cyl_indices); 

% Main loop
disp(['Simulation started, running ' num2str(timesteps) ' timesteps...']);
run_LBM_loop(f, u, v, rho, omega, u_lb, b_cyl_indices, i_cyl_indices, a_cyl_indices, boundary_links, ... 
             dh, dt, timesteps, calc_coeff, sample_factor, x_cyl, y_cyl, cyl_rad_nodes, x_sampled, y_sampled, ...
             update_every_iter, vis, show_vector_field);
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

%%%% DISPLAYING STREAMLINES [note: idk if streamlines is protected word] 
% Streamlines 
streamlines = streamline(x_sampled, y_sampled, sample_u, sample_v, xstart, ystart); 
set(streamlines, 'color', 'black'); 

disp('**************Done!****************');
disp('Done!');
%} 