function PotEng1 = ExactPotEnergy1(r, Poisson, E, rho, Omega, Ri, Ry)
A = ((3+Poisson)/8)*rho*(Omega^2)*(Ry^2 - Ri^2);
B = ((3+Poisson)/8)*rho*(Omega^2)*(Ri^2)*(Ry^2);
sigma_r = A - B./(r.^2) - ((3+Poisson)/8)*rho*Omega^2*r.^2;
sigma_phi = A + B./(r.^2) - ((1+3*Poisson)/8)*rho*Omega^2*r.^2;
Epsilon_r = (sigma_r - Poisson*sigma_phi)/(E*(1+Poisson));
Epsilon_phi = (sigma_phi - Poisson*sigma_r)/(E*(1+Poisson));
PotEng1 = (1/2) * (Epsilon_r .* sigma_r + Epsilon_phi .* sigma_phi) .* r;
end