function [R, U, Numsigmar, Numsigmaphi, K, Fext] = Numsol(E, rho, Ri, Ry, Omega, Poisson, N)

%Preallocate
K = zeros(N);
Fext = zeros(N,1);
r = linspace(Ri, Ry, N);
Numsigmar = zeros(N,1);
Numsigmaphi = zeros(N,1);
R = zeros(1,N);