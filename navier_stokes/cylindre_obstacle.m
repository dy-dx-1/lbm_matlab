% A Lattice Boltzmann (single relaxation time) D2Q9 solver,
% with the Spalart Allmaras turbulence model, on a lid-driven cavity. 
% Cell centers (nodes) are placed on the boundaries. 
% Author: Robert Lee
% Email: rlee32@gatech.edu

% Modifications du code effectuee dans le cadre du cours MEC3900-Projet Integrateur III
% faites par Nicolas Sarabia-Benoit, a Polytechnique Montreal 

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
nx = 135; % total # of nodes used in the sim, 16641 corresponds to 129x129 for GHIA comparison 
duct_ratio = 3; % ratio for the rectangel such that Length = ratio*Height 
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
y_cyl = round(ny/2)+1; % Y position of center of cyl, slightly offset 
cyl_matrix = generate_obstacle_matrix(X, Y, x_cyl, y_cyl, cyl_rad_nodes, 'circle');  % Matrix where 1 represents a cylinder node 
cyl_indices = find(cyl_matrix); % linear indexation of non zero elemetns, will be used to apply BB https://www.mathworks.com/help/matlab/ref/find.html

% prepping calculations for lift and drag coeff (see "obstacles/aero_coeffs.m")
calc_coeff = (2/(rho0*(u_lb^2)*2*cyl_rad_nodes)); 
% Displaying info 
display_sim_info(nx, ny, timesteps, Re, dh, dt, nu_lb, tau);

% Initialize.
rho = rho0*ones(ny,nx);
u = zeros(ny,nx);
v = zeros(ny,nx);
u(2:end-1,1) = u_lb;

f = compute_feq(rho,u,v);
% Apply meso BCs.
f = apply_meso_obs(f, u_lb, cyl_indices); 

% Main loop.
disp(['Running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps
    if (mod(iter,round(timesteps/10))==0)
        disp(['Running ... ' num2str(iter/timesteps*100) '% completed']);
    end
    
    % Collision.
    f = collide_mrt(f, u, v, rho, omega);
    
    % Apply meso BCs.
    f = apply_meso_obs(f, u_lb, cyl_indices); 

    % Streaming.
    f = stream(f);
    
    % Calculation of drag 
    [cd, cl] = aero_coeffs(f, cyl_indices, dh, dt, calc_coeff); 
    
    
    % Apply meso BCs.
    f = apply_meso_obs(f, u_lb, cyl_indices); 
    
    % Apply macro variables
    [u,v,rho] = apply_macro_obs(f, u_lb); 
    
    %{
    %AVERAGE DENSITY VISUALISATION
    average_density = mean(rho, 'all');
    fileID = fopen('C:/Users/Nicolas/Downloads/average_density.txt', 'a');
    fprintf(fileID, '%f\n', average_density); 
    fclose(fileID); 
    disp(average_density); % displaying average density to check for conservation
    %} 
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
        uu(cyl_indices) = nan; 
        imagesc(flipud(uu));
        % rectangle function is easiest to draw a circle, pos vector outlines lower left corner and height and width, curvature makes it a circle 
        rectangle('Position', [x_cyl-cyl_rad_nodes y_cyl-cyl_rad_nodes-1 cyl_rad_nodes*2 cyl_rad_nodes*2], 'Curvature', [1 1], 'FaceColor', 'red')
        colorbar
        axis equal; 
        
        %{
        subplot(2,1,2); 
        % Getting pressure field & displaying it 
        p = rho.*(1/3)*((dh/dt))^2;      
        contourf(p, 50);
        colorbar
        axis equal; 
        %} 
        drawnow
    %end
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