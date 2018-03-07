%% Calculating the potential energy
function [E_pot_num, E_pot_an] = E_sol(Ri, Ry, E, pois, Omega, rho, N_E_pot)
    E_pot_num = zeros(1,N_E_pot-1);

    for N = 2:N_E_pot
        [~, u_num, ~, ~, K, f_an] = Num_sol(E, rho, Ri, Ry, Omega, pois, N);
        E_pot_num(N-1) = 0.5.*u_num'*K*u_num - u_num'*f_an;
    end

    E_pot_an = integral(@(r) E_p_an1(r, pois, E, rho, Omega, Ri, Ry),Ri, Ry)...
                -integral(@(r) E_p_an2(r, rho, Omega, pois, E, Ri, Ry), Ri, Ry);
end

%% Internal functions
% Potential energy analytical 1
function E_p1 = E_p_an1(r, pois, E, rho, Omega, Ri, Ry)
    A = ((3+pois)/8)*rho*(Omega^2)*(Ry^2 - Ri^2);
    B = ((3+pois)/8)*rho*(Omega^2)*(Ri^2)*(Ry^2);
    sig_r = A - B./(r.^2) - ((3+pois)/8)*rho*Omega^2*r.^2;
    sig_phi = A + B./(r.^2) - ((1+3*pois)/8)*rho*Omega^2*r.^2;
    eps_r = (sig_r - pois*sig_phi)/(E*(1+pois));
    eps_phi = (sig_phi - pois*sig_r)/(E*(1+pois));
    E_p1 = (1/2) * (eps_r .* sig_r + eps_phi .* sig_phi) .* r;
end

% Potential energy analytical 2
function E_p2 = E_p_an2(r, rho, Omega, pois, E, Ri, Ry)
    u = ((3+pois)/8)*(1-pois)*((rho*Omega^2)/E).*r.*(Ri^2+Ry^2-((1+pois)/(3+pois)).*r.^2+((1+pois)/(1-pois))*((Ri^2*Ry^2)./r.^2));
    E_p2 = u .* rho * Omega^2 .* r.^2;
end