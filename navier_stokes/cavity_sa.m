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
nodes = 129;
timesteps = 1000000; 
nutilde0 = 1e-5; % initial nutilde value (should be non-zero for seeding).
u_lb = 0.1; 

%%% Simulation parameters.
Re = 3200; % Nombre de Reynolds, a  commenter pour imposer viscosite cinematique 
%nu_p = 1.2e-3; % 1.586e-5; % Viscosite cinematique, commenter pour imposer Reynolds 
if (nu_p~=0) % Dans ce cas, nu_p n'a ete update, donc il n'est pas commente et il faut evaluer Re avec sa valeur 
    disp("nu_p impose, Re calcule a partir de nouvelle valeur de nu_p.");
else 
    disp("Re impose, nu_p impose a partir de Re"); 
    nu_p = L_p*U_p / Re; 
end
disp(strcat("nu_p = ", num2str(nu_p)));


dx_p = L_p/(N-1); % Espacement physique noeuds 
Cl = dx_p; % Coeff convertion tel que dx_p = Cl*dx_lb 
dt_p = dx_p*(u_lb/U_p); 
Ct = dt_p % Coeff convertion tel que dt_p = Cl*dt_lb 
Cnu = (Cl^2)/Ct % Coeff tel que nu_p = Cnu * nu_lb 

nu_lb = nu_p/Cnu
tau = 3*nu_lb + 0.5;
omega = 1 / tau;

%% testing zone 
dt = dt_p 
dh = dx_p 

% Displaying info 
disp(['Reynolds number: ' num2str(Re)]);
disp(['Lattice viscosity: ' num2str(nu_lb)]);
disp(['Original relaxation time: ' num2str(tau)]);
disp(['Physical relaxation parameter: ' num2str(omega)]);
disp(['Lattice speed: ' num2str(u_lb)]);

% Info sur le setup numerique et checks de stabilite 
total_time = dt_p*timesteps; 
disp(strcat("Total real simulation time : ", num2str(total_time), " s")); 
ratio_relax_dt = tau/dt_p; 
if (ratio_relax_dt<1)
    disp(strcat("WARNING! tau/dt ratio lower than 1, ref stab checks for BGK unavailable. tau/dt = ", num2str(ratio_relax_dt)))
end 
cond_stab_bgk = sqrt(2/3)*(dh/dt_p); % juste pour avoir un idee par rapport a  BGK
disp(strcat("BGK stability condition (should be >U_p): ", num2str(cond_stab_bgk))); 
if (U_p>cond_stab_bgk)
    disp("WARNING!!! BGK stab condition not respected"); 
end 
if (tau>5 && tau<0.5) 
    disp("WARNING!! Relaxation time tau not stable"); 
end 

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
    if (mod(iter,timesteps/10)==0)
        disp(['Ran ' num2str(iter) ' iterations']);
        % extracting velocity data along the middle 
        u_center = flipud(extractRowOrColumn(u, 'col', round(nodes/2)))/u_lb; % u along the vertical line at center 
        v_center = flipud(extractRowOrColumn(v, 'row', round(nodes/2)))/u_lb; % v along the horizontal line at center 
        % displaying velocities in array form 
        out_center_speeds = zeros(nodes, 2); 
        out_center_speeds(:,1) = u_center; 
        out_center_speeds(:,2) = v_center; 
        disp(out_center_speeds); 
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