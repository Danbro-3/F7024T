close all
clear
clc

syms u E omega v rho f_r
constex1=1; %acts as placeholders
constex2=2; %

R_i=constex1;                       %inner radii
R_o=constex2;                       %outer radii
r=linspace(R_i,R_o);                %radius
v = constex1;                       %velocity
omega= constex1;                    %angular velocity, rad/s

%u=                     %Nodal value
%E=                     %elastic modulus
%v=                     %poissons ratio
%rho=                   %density

A = (3+v)/8*rho*omega^2*(Ri^2*Ro^2-Ri^4)/Ri^2;
B = (3+v)/8*rho*omega^2*Ri^2*Ro^2;

%eps_r = diff(u);
%eps_phi = u/r;
%S_r = (E/(1-v^2))*[1 v]*[ers_r;eps_phi]
%S_phi = (E/(1-v^2))*[v 1]*[ers_r;eps_phi]

%f_r=rho*r*omega^2      %volumetric loading for rotating disc
%L = @(u)diff(u,r,2)+diff(u,r)*./r-u/r^2+(1-v^2)*f_r/E;

%--------------------------------------------------
% Galerkin and Ritz


