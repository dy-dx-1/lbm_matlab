function f = inlet_bc(f, u_lb, side)
% D2Q9
% Applies an inlet BC to the distirbution functions in 2d matrix form.
% Assumes boundary nodes are located on the boundary (as opposed to half-cell away).

if strcmp(side, 'west') % West inlet, fixed velocity. BB method 
    %rho_boundary = 1; % densite constante, ecoulement incompressible 
    [~, ~, rho] = reconstruct_macro_all(f); 
    rho_boundary = 0.5*(rho(2:end-1, end) + rho(2:end-1, end-1)); 
    cs = sqrt(1/3); 

    % poids et directions
    w_4 = 1/9;
    w_8 = 1/36; 
    w_7 = 1/36;  
    
    c_4 = [-1; 0];
    c_8 = [-1; -1];
    c_7 = [-1; 1];

    % vitesse, on prend la vitesse imposée u_lb 
    [rows, ~, ~] = size(f); 
    % u_w size [nodes-2X2], each row has vector u,v. -2 because the populations are indexed 2:end-1
    u_w = repmat([u_lb, 0], rows-2, 1);

    f(2:end-1,1,2) = f(2:end-1,1,4) - ((2*w_4*rho_boundary/(cs^2)).*(u_w*c_4));
    f(2:end-1,1,6) = f(2:end-1,1,8) - ((2*w_8*rho_boundary/(cs^2)).*(u_w*c_8));
    f(2:end-1,1,9) = f(2:end-1,1,7) - ((2*w_7*rho_boundary/(cs^2)).*(u_w*c_7));
end 
%{
if strcmp(side, 'west') % West inlet, fixed velocity.    Zou/he 
    rho_west = ( 1 / ( 1 - u_lb ) ) * ...
        ( f(2:end-1,1,1) + f(2:end-1,1,3) + f(2:end-1,1,5) + ...
        2*( f(2:end-1,1,4) + f(2:end-1,1,7) + f(2:end-1,1,8) ) );
    f(2:end-1,1,2) = f(2:end-1,1,4) + 2 / 3 * u_lb * rho_west;
    f(2:end-1,1,6) = f(2:end-1,1,8) + u_lb / 6 * rho_west;
    f(2:end-1,1,9) = f(2:end-1,1,7) + u_lb / 6 * rho_west;
end
%}