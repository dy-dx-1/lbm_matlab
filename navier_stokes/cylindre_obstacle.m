% Code simulating turbulent flow around a cylinder using the Lattice Boltzmann Method. 
% This code is based on the code in "cavity_sa.m" by Robert Lee (rlee32@gatech.edu)
% and modified by Nicolas Sarabia-Benoit (nicolas.sarabia-benoit@polymtl.ca) in the context
% of the course MEC3900 at Polytechnique Montreal.

clear;close all;clc;fclose('all'); 
delete 'C:/Users/Nicolas/Desktop/PI3_calcs/results_coeffs/CD.txt'
delete 'C:/Users/Nicolas/Desktop/PI3_calcs/results_coeffs/CL.txt'
delete 'C:/Users/Nicolas/Desktop/PI3_calcs/results_coeffs/CP.txt'
addpath basic
addpath bc
addpath turbulence
addpath verif_assets
addpath obstacles

%%% Base parameters
Re = 900;         % Reynolds number
tau = 0.62;       % Relaxation time
total_time = 5;   % Total simulation time
rho0 = 1;         % Initial density (adim)

%%% Geometric parameters.
nx = 250;             % # of nodes in the x direction
duct_ratio = 1;       % ratio for the rectangular domain such that Length = ratio*Height 
cyl_size_ratio = 0.1; % Diam of cyl as a fraction of the duct height
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
x_cyl = round(nx/5);                          % X position of center of cyl 
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
pressure_calc_coeff = (1/3)*((dh/dt)^2); 
p_inf = rho0*pressure_calc_coeff; % Far upstream pressure (reference for Cp)
p_divider = 0.5*rho0*u_lb*u_lb; % Such that Cp = (p-p_inf)/p_divider 

%%% Prepping for visualization (to avoid loop calculations)
% Velocity vector field (quiver)
sample_factor = 8;             % To have a less dense quiver field
[x_sampled, y_sampled] = meshgrid(1:sample_factor:nx, 1:sample_factor:ny); 
% Streamlines
N = 75; % number of seed locations ; used to display streamlines 
xstart = nx*rand(N,1); 
ystart = ny*rand(N,1);

%%% Visualization choices 
% This is a cell array that contains the visualization options for the simulation
% It needs to respect the following format: {aero_coeffs, density, velocity||pressure||velocity_pressure, show_shape}
% To not call a function use @placeholder
% Functions are defined in /verif_assets/
vis = {@show_aero_coeffs, @placeholder, @show_velocity_and_pressure, @placeholder};
update_every_iter = 0; % 1 to update every iteration, 0 to update every 10% of the simulation, use 0 for large simulations
show_vector_field = 1; % 1 to show velocity quiver, 0 to not 

%%% Simulation ---------------------------------------------------------------------------------------------------------------------
% Initialization
rho = rho0*ones(ny,nx);
u = zeros(ny,nx);
v = zeros(ny,nx);
u(2:end-1,1) = u_lb;

f = compute_feq(rho,u,v);
f = apply_meso_obs(f, u_lb, b_cyl_indices, boundary_links); 

% Main loop
disp(['Simulation started, running ' num2str(timesteps) ' timesteps...']);
[u, v, rho, f] = run_LBM_loop(f, u, v, rho, omega, u_lb, b_cyl_indices, i_cyl_indices, a_cyl_indices, boundary_links, ... 
             dh, dt, pressure_calc_coeff, p_inf, p_divider, timesteps, calc_coeff, sample_factor, x_cyl, y_cyl, cyl_rad_nodes, x_sampled, y_sampled, ...
             update_every_iter, vis, show_vector_field);

%%% Post simulation calculations ----------------------------------------------------------------------------------------------
% Displaying speeds on the centerlines 
%display_center_speeds(u, v, u_lb, nx, ny); 

% Last visual update of iters 
vis{3}(show_vector_field, u, v, u_lb, sample_factor, x_sampled, y_sampled, rho, pressure_calc_coeff, a_cyl_indices);
fig = findobj('type', 'figure'); % Find all figure objects
if numel(fig) == 2 % Check if there are two figures
    ax1 = fig(1).Children; % Axes of the first figure
    ax2 = fig(2).Children; % Axes of the second figure

    % Formatting for the first figure
    title(ax1, "Champ de pression");
    xlabel(ax1, "Noeuds en X");
    ylabel(ax1, "Noeuds en Y");
    colorbar(ax1);
    axis(ax1, 'equal');

    % Formatting for the second figure
    title(ax2, "Champ de vitesse");
    xlabel(ax2, "Noeuds en X");
    ylabel(ax2, "Noeuds en Y");
    colorbar(ax2);
    axis(ax2, 'equal');
else
    xlabel("Noeuds en X");
    ylabel("Noeuds en Y");
    colorbar;
    axis equal;
end
% Streamlines 
%show_streamlines(x_sampled, y_sampled, u, v, sample_factor, xstart, ystart); 

disp('**************Done!****************');