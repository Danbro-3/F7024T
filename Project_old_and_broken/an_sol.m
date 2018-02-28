%% Analytical solution of disc problem
function [r, u, sig_r, sig_phi] = an_sol(E, rho, Ri, Ry, Omega, pois, N)
% Analytical solution of disc problem
    r = linspace(Ri, Ry, N);
    u = zeros(1, N);
    sig_r = zeros(1, N);
    sig_phi = zeros(1, N);
    A = ((3+pois)/8)*rho*(Omega^2)*(Ry^2 - Ri^2);
    B = ((3+pois)/8)*rho*(Omega^2)*(Ri^2)*(Ry^2);
    for i = 1:N
        u(i) = ((3+pois)/8)*(1-pois)*((rho*Omega^2)/E)*r(i)*(Ri^2 + Ry^2-((1+pois)/(3+pois))*r(i)^2+((1+pois)/(1-pois))*((Ri^2 * Ry^2)/r(i)^2));
        sig_r(i) = A-B/(r(i)^2)-((3+pois)/8)*rho*Omega^2*r(i)^2;
        sig_phi(i) = A + B/(r(i)^2)-((1+3*pois)/8)*rho*Omega^2*r(i)^2;
    end
end