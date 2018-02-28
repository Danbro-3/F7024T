%% Setup 
% Clear old workspace
close all
clear, clc

Omega = 2*pi;   % rad/s
Ri = 0.1;       % m
Ry = 2.0;       % m
N = 500;        % nodes for analytical and numerical solution
N_E_pot = 50;   % nodes for potential energy

%% Material values: Platinum
E = 168E9;      % 015376 %[Pa] Elasticity, Youngs module
rho = 21450;    % [kg/m^3] Density
pois = 0.39; % Poissons ratio

%% Calculations
% Find the analytical solution, numerical solution and potential energy
[r_an, u_an, sig_r_an, sig_phi_an]          = an_sol(E, rho, Ri, Ry, Omega,pois, N);
[r_num, u_num, sig_r_num, sig_phi_num, K, ~]= Num_sol(E, rho, Ri, Ry, Omega,pois, N);
[E_pot_num, ~]                              = E_sol(Ri, Ry, E, pois, Omega, rho, N_E_pot);

%% Plotting
subplot(2,2,1)
plot(r_an, sig_r_an, r_num, sig_r_num, 'Linewidth', 1.5); 
xlabel('\bf{radius, $r$ [m]}', 'Interpreter','latex');
ylabel('\bf{stress, $\sigma_r$ [Pa]}', 'Interpreter','latex'); 
legend('Analytic method', 'Numerical method');
title('\bf{Stress in `$r$'' direction}', 'Interpreter','latex'); grid on
subplot(2,2,2) 
plot(r_an, sig_phi_an, r_num, sig_phi_num); 
xlabel('\bf{radius, $r$ [m]}', 'Interpreter','latex');
ylabel('\bf{stress, $\sigma_{\phi}$ [Pa]}', 'Interpreter','latex'); 
legend('Analytic method', 'Numerical method');
title('\bf{Stress in $\phi$ direction}', 'Interpreter','latex'); grid on
subplot(2,2,3)
plot(r_an, u_an, r_num, u_num); 
xlabel('\bf{radius, $r$ [m]}', 'Interpreter','latex');
ylabel('\bf{Displacement [m]}', 'Interpreter','latex'); 
legend('Analytic method', 'Numerical method', 'location', 'NorthWest');
title('\bf{Displacement over `$r$'' direction}', 'Interpreter','latex'); grid on
subplot(2,2,4)
plot(2:N_E_pot, E_pot_num);
xlabel('\bf{Number of nodes, n}', 'Interpreter','latex');
ylabel('\bf{energy, $E_P$}', 'Interpreter','latex');
title('\bf{Total potential energy (numeric)}', 'Interpreter','latex'); grid on