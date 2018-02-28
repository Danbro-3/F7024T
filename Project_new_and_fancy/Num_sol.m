%% Numerical solution of disc problem
function [R, U, sig_r_num, sig_phi_num, K, F_an] = Num_sol(E, rho, Ri, Ry, Omega, pois, N)
% Numerical solution of disc problem
    %% setup
    K = zeros(N);
    F_an = zeros(N,1);
    r = linspace(Ri, Ry, N);
    sig_r_num = zeros(N,1);
    sig_phi_num = zeros(N,1);
    R = zeros(1,N);
    %% formulas from compendium
    s = [-1/sqrt(3), 1/sqrt(3)];                        % From table (7.1)
    w = [1, 1];                                         % From table (7.1)
    Emod = [E/(1-pois^2), E*pois/(1-pois^2);
            E*pois/(1-pois^2), E/(1-pois^2)];           % From equation (8.6)    
        
    %% calc
    for i = 1:N-1
        k = zeros(1,2); f_an = zeros(2,1);
        B = [-1/(r(i+1)-r(i)), 1/(r(i+1)-r(i)); 0, 0];  % From equation (6.57)
        for id = 1:2 % node id
            [k, f_an] = f_an_comp(s(id), r(i), r(i+1), B, Emod, w(id), k, f_an, Omega, rho);
        end
        K(i,i) = K(i,i)+k(1,1);
        K(i,i+1) = K(i,i+1)+k(1,2);
        K(i+1,i) = K(i+1,i)+k(2,1);
        K(i+1,i+1) = K(i+1,i+1)+k(2,2);

        F_an(i) = F_an(i)+f_an(1);
        F_an(i+1) = F_an(i+1)+f_an(2);
    end

    U = K\F_an;

    for i = 1:N-1
        r_sig = zeros(1,2); sig_r = zeros(1,2); sig_phi = zeros(1,2);
        B = [-1/(r(i+1)-r(i)), 1/(r(i+1)-r(i)); 0, 0];  % From equation (6.57)
        u = [U(i); U(i+1)];
        for id = 1:2 % node id
            [sig_r(id), sig_phi(id), r_sig(id), B] = sig_comp_num(s(id), r(i), r(i+1), B, Emod, u);
        end
        
        R(i) = r_sig(1);
        R(i+1) = r_sig(2);
        sig_r_num(i) = sig_r(1);
        sig_r_num(i+1) = sig_r(2);
        sig_phi_num(i) = sig_phi(1);
        sig_phi_num(i+1) = sig_phi(2);
    end
end

%% Internal functions
%Computating k and f_an stuff
function [k, f_an] = f_an_comp(s, r1, r2, B, E, w, k, f_an, Omega, rho)
    N = [(1-s)/2 , (1+s)/2];
    r = r1*N(1) + r2*N(2);
    B(2,:) = [N(1)/r ; N(2)/r];
    k = k + B'*E*B*r*w;
    f_an = f_an + N'*rho*(Omega^2)*(r^2)*w;
end

%Computating sigma stuff
function [sig_r, sig_phi, r, B] = sig_comp_num(s, r1, r2, B, E, u)
    N = [(1-s)/2, (1+s)/2];
    r = N(1)*r1+N(2)*r2;
    B(2,:) = [N(1)/r , N(2)/r];
    epsilon = B*u;
    sig = E*epsilon;
    sig_r = sig(1);
    sig_phi = sig(2);
end
