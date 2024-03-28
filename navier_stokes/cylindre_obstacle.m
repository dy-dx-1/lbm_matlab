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
Re = 100; 
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
y_cyl = round(ny/2); % Y position of center of cyl
cyl_matrix = generate_obstacle_matrix(X, Y, x_cyl, y_cyl, cyl_rad_nodes, 'circle');  % Matrix where 1 represents a cylinder node 
boundary_cyl_matrix = mark_boundary_nodes(cyl_matrix); % this code uses halfway BB, so the boundary nodes are technically still fluid nodes
boundary_links = find_boundary_links(cyl_matrix-boundary_cyl_matrix); % need to use only the 'inside without boundary' because halfwayBB makes it so outher boundary is 'fluid' & in find_boundary_links we use the 0's to find links
% linear indexation of non zero elemetns, will be used to apply BB https://www.mathworks.com/help/matlab/ref/find.htm
a_cyl_indices = find(cyl_matrix); % boundary and inside 
b_cyl_indices = find(boundary_cyl_matrix);  % only boundary
i_cyl_indices = find(cyl_matrix-boundary_cyl_matrix);  % only inside without boundary

% Displaying simulation info 
display_sim_info(nx, ny, timesteps, Re, dh, dt, nu_lb, tau);

% prepping calculations for lift and drag coeff (see "obstacles/aero_coeffs.m")
% and data needed for visualization 
calc_coeff = (2/(rho0*(u_lb^2)*2*cyl_rad_nodes)); 
sample_factor = 3; % To have a less dense quiver field
[x_sampled, y_sampled] = meshgrid(1:sample_factor:nx, 1:sample_factor:ny);
x_sampled = flipud(x_sampled); 
y_sampled = flipud(y_sampled); 

N = 75; % number of seed locations
xstart = nx*rand(N,1); 
ystart = ny*rand(N,1);

% Initialize.
rho = rho0*ones(ny,nx);
u = zeros(ny,nx);
v = zeros(ny,nx);
u(2:end-1,1) = u_lb;

f = compute_feq(rho,u,v);
f = apply_meso_obs(f, u_lb, b_cyl_indices); 

% Main loop.
disp(['Simulation started, running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps    
    % Collision.
    f = collide_mrt(f, u, v, rho, omega);
    
    % Apply meso BCs.
    f = apply_meso_obs(f, u_lb, b_cyl_indices); 

    % Calculation of drag and lift coefficient 
    [cd, cl] = aero_coeffs(f, b_cyl_indices, boundary_links, dh, dt, calc_coeff); 
    disp(cd); 
    
    % Streaming.
    f = stream(f);    
    
    % Apply meso BCs.
    f = apply_meso_obs(f, u_lb, b_cyl_indices); 
    
    % Apply macro variables
    [u,v,rho] = apply_macro_obs(f, u_lb, i_cyl_indices); 
    
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
    
    % VISUALIZATION & progress tracking 
    if (mod(iter,round(timesteps/10))==0)
        disp(['Running ... ' num2str(iter/timesteps*100) '% completed']);
    end
        %subplot(2,1,1); 
        uu = sqrt(u.^2+v.^2) / u_lb;
        uu(a_cyl_indices) = nan; 
        imagesc(flipud(uu));
        % Arrow vector field 
        hold on;
        sample_u = flipud(u(1:sample_factor:end, 1:sample_factor:end));
        sample_v = flipud(v(1:sample_factor:end, 1:sample_factor:end));
        quiver(x_sampled, y_sampled, sample_u, sample_v, 'r'); 
        hold off; 

        % rectangle function is easiest to draw a circle, pos vector outlines lower left corner and height and width, curvature makes it a circle 
        %rectangle('Position', [x_cyl-cyl_rad_nodes y_cyl-cyl_rad_nodes-1 cyl_rad_nodes*2 cyl_rad_nodes*2], 'Curvature', [1 1], 'FaceColor', 'red')
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
    %end UNCOMMENT AND COMMENT THE NEXT ONE UP FOR LIMITED VISUALISATION
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

%%%% DISPLAYING STREAMLINES [note: idk if streamlines is protected word] 
% Streamlines 
streamlines = streamline(x_sampled, y_sampled, sample_u, sample_v, xstart, ystart); 
set(streamlines, 'color', 'black'); 

disp('**************Done!****************');
disp('Done!');
%} 