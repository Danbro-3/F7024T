function [k, fext] = kfextComp(s, r1, r2, B, E, w, k, fext, Omega, rho)
N = [(1-s)/2 , (1+s)/2];
r = r1*N(1) + r2*N(2);
B(2,:) = [N(1)/r ; N(2)/r];
k = k + B'*E*B*r*w;
fext = fext + N'*rho*(Omega^2)*(r^2)*w;
end