function [sigmar, sigmaphi, r] = sigmaComp(r1, r2, u, E, Poisson)
s = [-1/sqrt(3), 1/sqrt(3)]; % p 49 in comp
r = zeros(1,2);
sigmar = zeros(1,2);
sigmaphi = zeros(1,2);
Emod = [E/(1-Poisson^2), E*Poisson/(1-Poisson^2); E*Poisson/(1-Poisson^2), E/(1-Poisson^2)]; % E = sigma/epsilon, from eq. (8.6)
B = [-1/(r2-r1), 1/(r2-r1) ;0 , 0]; % From eq. (6.57)
for id = 1:2 % two nodes
[sigmar(id), sigmaphi(id), r(id), B] = NumsigmaComp(s(id), r1, r2, B, Emod, u);
end
end