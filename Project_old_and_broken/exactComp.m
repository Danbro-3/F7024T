function [u, sigmar, sigmaphi] = exactComp(r, E, Poisson, rho, Omega, Ri, Ry, A, B)
u = ((3+Poisson)/8)*(1-Poisson)*((rho*Omega^2)/E)*r*(Ri^2 + Ry^2-((1+Poisson)/(3+Poisson))*r^2+((1+Poisson)/(1-Poisson))*((Ri^2 * Ry^2)/r^2));
sigmar = A-B/(r^2)-((3+Poisson)/8)*rho*Omega^2*r^2;
sigmaphi = A + B/(r^2)-((1+3*Poisson)/8)*rho*Omega^2*r^2;
end