% Lid-driven cavity with a cut corner.
% A Lattice Boltzmann D2Q9 solver.
% This features a non-lattice-aligned wall! 
% Cell centers (nodes) are placed on the boundaries. 
% Author: Robert Lee
% Email: rlee32@gatech.edu

clear;close all;clc;

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
L_p = 0.6;%1.1; % Cavity dimension. 
U_p = 6;%1.1; % Cavity lid velocity.
nu_p = 1.2e-3;%1.586e-5; % Physical kinematic viscosity.
rho0 = 1;
cut_start_y = 0.5; % non-dimensional y-position on the west boundary.
cut_end_x = 0.5; % non-dimensional x-position on the south boundary.
% Discrete/numerical parameters.
nodes = 100;
dt = .002;
timesteps = 10000;

% Derived nondimensional parameters.
Re = L_p * U_p / nu_p;
disp(['Reynolds number: ' num2str(Re)]);
% Derived physical parameters.
t_p = L_p / U_p;
disp(['Physical time scale: ' num2str(t_p) ' s']);
% Derived discrete parameters.
dh = 1/(nodes-1);
nu_lb = dt / dh^2 / Re;
disp(['Lattice viscosity: ' num2str(nu_lb)]);
tau = 3*nu_lb + 0.5;
disp(['Relaxation time: ' num2str(tau)]);
omega = 1 / tau;
disp(['Relaxation parameter: ' num2str(omega)]);
u_lb = dt / dh;
disp(['Lattice speed: ' num2str(u_lb)])

% Determine which lattice vectors are relevant to the cut.
parallel = [-cut_end_x, cut_start_y];
cut_length = norm(parallel);
unit_parallel = parallel / cut_length;
unit_normal = [-parallel(1), parallel(2)] / cut_length;
pgram_height = cut_length * dt;
c = zeros(9,2);
c(1,:) = [0, 0];
c(2,:) = [1, 0];
c(3,:) = [0, 1];
c(4,:) = [-1, 0];
c(5,:) = [0, -1];
c(6,:) = [1, 1];
c(7,:) = [-1, 1];
c(8,:) = [-1, -1];
c(9,:) = [1, -1];
valid = zeros(9,1);
for k = 1:9
    valid(k) = dot(unit_normal,c(k,:)) < 0;
end
c_wall = zeros(sum(valid),2);
counter = 1;
for k = 1:9
    if valid(k)
        c_wall(counter, :) = c(k, :);
        counter = counter + 1;
    end
end
% Pgram defined by pgram_height, cut_length, unit_normal, unit_parallel.
touched = zeros(nodes,nodes,1);
coord_min = dh*(cumsum(ones(nodes,1))-1) - dh/2;
for j = 1:nodes
    for i = 1:nodes
        
    end
end

% Initialize.
f = ones(nodes,nodes,9);
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

% Main loop.
disp(['Running ' num2str(timesteps) ' timesteps...']);
for iter = 1:timesteps
    if (mod(iter,timesteps/10)==0)
        disp(['Ran ' num2str(iter) ' iterations']);
    end
    
    % Collision.
    f = collide(f, u, v, rho, omega);
    
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
    
    % VISUALIZATION
    % Modified from Jonas Latt's cavity code on the Palabos website.
    if (mod(iter,10)==0)
        uu = sqrt(u.^2+v.^2) / u_lb;
        imagesc(flipud(uu));
        colorbar
        axis equal off; drawnow
    end
end
disp('Done!');



