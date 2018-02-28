function [sigmar, sigmaphi, r, B] = NumsigmaComp(s, r1, r2, B, E, u)
N = [(1-s)/2, (1+s)/2];
r = N(1)*r1+N(2)*r2;
B(2,:) = [N(1)/r , N(2)/r];
epsilon = B*u;
sigma = E*epsilon;
sigmar = sigma(1);
sigmaphi = sigma(2);
end