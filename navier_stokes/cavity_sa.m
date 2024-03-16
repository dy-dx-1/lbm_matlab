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
L_p = 4;   % Cavity dimension. 
nu_p = 0.1; % Viscosite cinematique 
rho0 = 1;  % Densite initiale 

nodes = 129;
total_time = 20; % temps de simulation total en sec 
nutilde0 = 1e-5; % initial nutilde value (should be non-zero for seeding).

%%% Simulation parameters.
Re = 400; % Nombre de Reynolds, a  commenter pour imposer viscosite cinematique 

dx_p = L_p/(nodes-1); 
% dt_p chosen following diffusive scaling (p.278 Krüger book) such as dt_p proportional dx_p^2 
%alpha = u_lb/dx_p;  % u_lb = dt/dx --> u_lb = alpha*dx
%dt_p = alpha*dx_p*dx_p;

%%% testing
tau = 0.809;
u_lb = 0.03; 
nu_lb = (tau-0.5)/3;  % 0.103 avec tau = 0.809
dt_p = (nu_p/nu_lb)*u_lb*u_lb; 

U_p = Re*nu_p/L_p; %Vitesse de la paroi en m/s 
omega = 1/tau;

timesteps = round(total_time/dt_p); 

% Setting calculation params
dt = dt_p; 
dh = dx_p; 

% Displaying info 
display_sim_info(L_p, U_p, nodes, timesteps, Re, dh, dt, u_lb, tau);

% Determine macro variables and apply macro BCs
% Initialize macro, then meso.
rho = rho0*ones(nodes,nodes);
u = zeros(nodes,nodes);
v = zeros(nodes,nodes);
u(end,2:end-1) = u_lb;
% Initialize.
f = compute_feq(rho,u,v);
% Apply meso BCs.
f = moving_wall_bc(f,'north',u_lb);
f = wall_bc(f,'south');
f = wall_bc(f,'east');
f = wall_bc(f,'west');
% Initialize turbulence stuff.
d = compute_wall_distances(nodes);
nutilde = nutilde0*ones(nodes,nodes);
[omega, nut, nutilde] = update_nut(nutilde,nu_lb,dt,dh,d,u,v);

% Main loop.
disp(['Running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps
    if (mod(iter,round(timesteps/10))==0)
        disp(['Ran ' num2str(iter) ' iterations out of ' num2str(timesteps) '---']);
    end
    
    % Collision.
    f = collide_sa(f, u, v, rho, omega);
    
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
    [omega, nut, nutilde] = update_nut(nutilde,nu_lb,dt,dh,d,u,v);
    
    % Sanity check that the simulation is working, i.e checking that all non boundary nodes in the u matrix are not NaN
    if (any(isnan(u(2:end-1,2:end-1))))
        disp('Error: NaNs in u matrix. Exiting...');
        return;
    end

    % VISUALIZATION
    % Modified from Jonas Latt's cavity code on the Palabos website.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%COMMENT/UNCOMMENT TO SPEED UP SIMULATION
    %{
    if (mod(iter,10==0))
        uu = sqrt(u.^2+v.^2) / u_lb;
        imagesc(flipud(uu));
%       imagesc(flipud(nut)); %%% turbulence viscosity 
%       imagesc(flipud(omega));

        colorbar
        axis equal; 
        drawnow
    end
    %}
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