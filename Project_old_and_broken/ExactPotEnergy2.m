function PotEng2 = ExactPotEnergy2(r, rho, Omega, Poisson, E, Ri, Ry)
u = ((3+Poisson)/8)*(1-Poisson)*((rho*Omega^2)/E).*r.*(Ri^2+Ry^2-((1+Poisson)/(3+Poisson)).*r.^2+((1+Poisson)/(1-Poisson))*((Ri^2*Ry^2)./r.^2));
PotEng2 = u .* rho * Omega^2 .* r.^2;
end