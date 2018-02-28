function [R, U, Numsigmar, Numsigmaphi, K, Fext] = Numsol(E, rho, Ri, Ry, Omega,Poisson, N)
K = zeros(N);
Fext = zeros(N,1);
r = linspace(Ri, Ry, N);
Numsigmar = zeros(N,1);
Numsigmaphi = zeros(N,1);
R = zeros(1,N);
for Nelem = 1:N-1
[k, fext] = FEMComp(r(Nelem), r(Nelem+1), E, Poisson, Omega, rho);
K(Nelem,Nelem) = K(Nelem,Nelem)+k(1,1);
K(Nelem,Nelem+1) = K(Nelem,Nelem+1)+k(1,2);
K(Nelem+1,Nelem) = K(Nelem+1,Nelem)+k(2,1);
K(Nelem+1,Nelem+1) = K(Nelem+1,Nelem+1)+k(2,2);
Fext(Nelem) = Fext(Nelem)+fext(1);
Fext(Nelem+1) = Fext(Nelem+1)+fext(2);
end
U = K\Fext;
for Nelem = 1:N-1
u = [U(Nelem); U(Nelem+1)];
[sigmar, sigmaphi, rsigma] = sigmaComp(r(Nelem), r(Nelem+1), u, E, Poisson);
R(Nelem) = rsigma(1);
R(Nelem+1) = rsigma(2);
Numsigmar(Nelem) = sigmar(1);
Numsigmar(Nelem+1) = sigmar(2);
Numsigmaphi(Nelem) = sigmaphi(1);
Numsigmaphi(Nelem+1) = sigmaphi(2);
end
end