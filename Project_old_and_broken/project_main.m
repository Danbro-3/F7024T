
% Material values: Aluminium
E = 68.950E9;       % 015376 %[Pa] Elasticity, Youngs module
rho = 2.7E3;        % [kg/m^3] Density
Poisson = 0.334;    % Poissons ratio

Omega = 2*pi; % rad/s
Ri = 0.1; % m
Ry = 2.0; % m
N = 500;
N_exact = 100;
N_Errornodes = 18;
N_PotEnergy = 50;

[r_exact, u_exact, sigmar_exact, sigmaphi_exact] = exactsol(E, rho, Ri, Ry, Omega,Poisson, N_exact);
[r_num, u_num, sigmar_num, sigmaphi_num, K, Fext] = Numsol(E, rho, Ri, Ry, Omega,Poisson, N);

figure; 
plot(r_exact, u_exact, r_num, u_num); 
xlabel('r [m]'); 
ylabel('Displacement[m]'); 
legend('Analytic method', 'Numerical method');

figure; 
plot(r_exact, sigmar_exact, r_num, sigmar_num); 
xlabel('r [m]');
ylabel('\sigma_r [Pa]'); 
legend('Analytic method', 'Numerical method');
figure; 
plot(r_exact, sigmaphi_exact, r_num, sigmaphi_num); 
xlabel('r [m]');
ylabel('\sigma_{\phi} [Pa]'); 
legend('Analytic method', 'Numerical method');

[PotEnergy_num, PotEnergy_exact] = Energysol(Ri, Ry, E, Poisson, Omega, rho, N_PotEnergy);
