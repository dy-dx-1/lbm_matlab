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

%%% Physical base parameters.
L_p = 1;   % Cavity dimension. 
U_p = 1;   % Cavity lid velocity.
rho0 = 1;  % Densite initiale 
nu_p = 0;  % sert a  setup le choix d'imposer viscosite ou nombre de Reynolds 

nodes = 16641; % total # of nodes used in the sim, 16641 corresponds to 129x129 for GHIA comparison 
duct_ratio = 3; % ratio for the rectangel such that Length = ratio*Height

total_time = 20; %temps de simulation total en sec 
nutilde0 = 1e-5; % initial nutilde value (should be non-zero for seeding).
u_lb = 0.18; 

%%% Simulation parameters.
Re = 3200; % Nombre de Reynolds, a  commenter pour imposer viscosite cinematique 
%nu_p = 1.2e-3; % 1.586e-5; % Viscosite cinematique, commenter pour imposer Reynolds 
if (nu_p~=0) % Dans ce cas, nu_p n'a pas ete update, donc il n'est pas commente et il faut evaluer Re avec sa valeur 
    disp("nu_p impose, Re calcule a partir de nouvelle valeur de nu_p.");
else 
    disp("Re impose, nu_p impose a partir de Re"); 
    nu_p = L_p*U_p / Re; 
end
disp(strcat("nu_p = ", num2str(nu_p)));

% Derived domain nodes ; considering dx=dy=dh
nx = round(((1-duct_ratio)+sqrt((duct_ratio-1)^2 + 4*duct_ratio*nodes))*0.5); 
ny = round(nodes/nx); 

% Lbm params s
dx_p = L_p/(nx-1); 
dt_p = u_lb*dx_p;
nu_lb = nu_p*dt_p/(dx_p*dx_p);
tau = 3*nu_lb + 0.5; 
omega = 1/tau;

timesteps = round(total_time/dt_p); 

%% testing zone 
dt = dt_p; 
dh = dx_p; 

% Displaying info 
display_sim_info(L_p, U_p, nodes, nx, ny, timesteps, Re, dh, dt, nu_lb, tau);

% Determine macro variables and apply macro BCs
% Initialize macro, then meso.
rho = rho0*ones(ny,nx);
u = zeros(ny,nx);
v = zeros(ny,nx);
u(2:end-1,1) = u_lb;

% Initialize.
f = compute_feq(rho,u,v);
% Apply meso BCs.
f = wall_bc(f,'north');
f = wall_bc(f,'south');
f = pressure_bc(f,'east');
f = inlet_bc(f, u_lb, 'west');
% Initialize turbulence stuff.
d = compute_wall_distances(ny, nx);
nutilde = nutilde0*ones(ny, nx);
[omega, nut, nutilde] = update_nut(nutilde,nu_lb,dt,dh,d,u,v);

% Main loop.
disp(['Running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps
    if (mod(iter,timesteps/10)==0)
        disp(['Ran ' num2str(iter) ' iterations']);
        % extracting velocity data along the middle 
        % TO CHECK AFTER RECTANGLE IMPLEMENTATION!!! 
        u_center = flipud(extractRowOrColumn(u, 'col', round(nx/2)))/u_lb; % u along the vertical line at center 
        v_center = flipud(extractRowOrColumn(v, 'row', round(ny/2)))/u_lb; % v along the horizontal line at center 
        % displaying velocities in array form 
        disp([u_center;v_center]); 
    end
    
    % Collision.
    f = collide_sa(f, u, v, rho, omega);
    
    % Apply meso BCs.
    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = pressure_bc(f,'east'); 
    f = inlet_bc(f, u_lb, 'west'); % constant entry speed at the left 

    % Streaming.
    f = stream(f);
    
    % Apply meso BCs.
    f = wall_bc(f,'north');
    f = wall_bc(f,'south');
    f = pressure_bc(f,'east');  
    f = inlet_bc(f, u_lb, 'west'); % constant entry speed at the left 
    
    % Determine macro variables
    [u,v,rho] = reconstruct_macro_all(f);
    % Macro BCs 
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

    [omega, nut, nutilde] = update_nut(nutilde,nu_lb,dt,dh,d,u,v);
    
    % VISUALIZATION
    % Modified from Jonas Latt's cavity code on the Palabos website.
    if (mod(iter,10==0))
        uu = sqrt(u.^2+v.^2) / u_lb;
        imagesc(flipud(uu));
%       imagesc(flipud(nut)); %%% turbulence viscosity 
%       imagesc(flipud(omega));

        colorbar
        axis equal; 
        drawnow
    end
end
disp('Done!');