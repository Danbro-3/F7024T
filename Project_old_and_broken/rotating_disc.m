function [r_ans, u_ans, sigmar_ans, sigmaphi_ans]=rotating_disc(n)
    %% Initiation
    close all
    % Material values: Aluminium
    E = 68.950E9;%015376 %[Pa] Elasticity, Youngs module
    rho = 2.7E3; % [kg/m^3] Density
    pois = 0.334; % Poissons ratio

    % Some unknown values
    Omega = 2*pi; % rad/s
    Ri = 0.01; % m
    Ry = 1.0; % m
    %n = 15;
    r_span = linspace(Ri, Ry, n);
    
    %% Solver function
    [r_ans, u_ans, sigmar_ans, sigmaphi_ans] = rot_solver(E, rho, r_span, Omega, pois);
    
    %% Plot result
    plot(r_ans, u_ans, 'k', 'Linewidth', 1.5);
    grid on
    title('\bf{Displacement}', 'Interpreter','latex')
    xlabel('\bf{Displacement, u, [m]}', 'Interpreter','latex');
    ylabel('\bf{r [m]}', 'Interpreter','latex');
    
    figure(2)
    plot(r_ans, sigmar_ans, 'k', 'Linewidth', 1.5);
    grid on
    title('\bf{Stress in r-direction}', 'Interpreter','latex')
    xlabel('\bf{r [m]}', 'Interpreter','latex');
    ylabel('\bf{Stress $\sigma_r$}', 'Interpreter','latex');
    
    figure(3)
    plot(r_ans, sigmaphi_ans, 'k', 'Linewidth', 1.5);
    grid on
    title('\bf{Stress in rotational direction}', 'Interpreter','latex')
    xlabel('\bf{r [m]}', 'Interpreter','latex');
    ylabel('\bf{Stress $\sigma_{\phi}$}', 'Interpreter','latex');
end

function [r, u, sigmar, sigmaphi] = rot_solver(E, rho, r, Omega, pois)
    %% Initiate variables
    u = zeros(1, length(r));
    sigmar = zeros(1, length(r));
    sigmaphi = zeros(1, length(r));
    
    %% Using boundary conditions
    A = ((3+pois)/8)*rho*(Omega^2)*(r(end)^2 - r(1)^2);
    B = ((3+pois)/8)*rho*(Omega^2)*(r(1)^2)*(r(end)^2);

    %% Get results
    for i = 1:length(r)
        u(i) = ((3+pois)/8)*(1-pois)*((rho*Omega^2)/E)*r(i)*(r(1)^2 + r(end)^2 - ((1-pois)/(3+pois))*r(i)^2 + ((1+pois)/(1-pois))*((r(1)^2 * r(end)^2)/r(i)^2));
        sigmar(i) = A - B/(r(i)^2) - ((3+pois)/8)*rho*Omega^2*r(i)^2;   
        sigmaphi(i) = A + B/(r(i)^2) - ((1+3*pois)/8)*rho*Omega^2*r(i)^2;
    end
end