function [r, u, sigmar, sigmaphi] = exactsol(E, rho, Ri, Ry, Omega, Poisson, N)
r = linspace(Ri, Ry, N);
u = zeros(1, N);
sigmar = zeros(1, N);
sigmaphi = zeros(1, N);
A = ((3+Poisson)/8)*rho*(Omega^2)*(Ry^2 - Ri^2);
B = ((3+Poisson)/8)*rho*(Omega^2)*(Ri^2)*(Ry^2);
for i = 1:N
    [u(i), sigmar(i), sigmaphi(i)] = exactComp(r(i), E, Poisson, rho, Omega, Ri, Ry, A, B);
end