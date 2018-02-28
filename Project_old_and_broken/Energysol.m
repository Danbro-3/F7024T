function [PotEnergy_num, PotEnergy_exact] = Energysol(Ri, Ry, E, Poisson, Omega, rho, N_PotEnergy)
PotEnergy_num = zeros(1,N_PotEnergy-1);
for N = 2:N_PotEnergy
[~, u_num, ~, ~, K, Fext] = Numsol(E, rho, Ri, Ry, Omega, Poisson, N);
PotEnergy_num(N-1) = 0.5.*u_num'*K*u_num - u_num'*Fext;
end
figure; plot(2:N_PotEnergy, PotEnergy_num); xlabel('Number of nodes'); ylabel('Total potential energy in the numeric case');
PotEnergy_exact = integral(@(r) ExactPotEnergy1(r, Poisson, E, rho, Omega, Ri, Ry),Ri, Ry)-integral(@(r) ExactPotEnergy2(r, rho, Omega, Poisson, E, Ri, Ry), Ri, Ry);
end