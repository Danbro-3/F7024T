function [k, fext] = FEMComp(r1, r2, E, Poisson, Omega, rho)
k = zeros(2);
fext = zeros(2,1);
s = [-1/sqrt(3), 1/sqrt(3)]; % p. 49 in comp
w = [1, 1]; % p. 49 in comp
Emod = [E/(1-Poisson^2), E*Poisson/(1-Poisson^2); E*Poisson/(1-Poisson^2), E/(1-Poisson^2)]; % E = sigma/epsilon, uttryck frn eq. (8.6)
B = [-1/(r2-r1), 1/(r2-r1); 0, 0]; % From eq. (6.57)
for id = 1:2 % Two nodes
[k, fext] = kfextComp(s(id), r1, r2, B, Emod, w(id), k, fext, Omega, rho);
end
end