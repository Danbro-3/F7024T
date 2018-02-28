function [ErrorTable] = WrongError(N_Errornodes, E, rho, Ri, Ry, Omega, Poisson)
format long
Feeeeel = zeros(1,N_Errornodes-1);
Cond = zeros(1,N_Errornodes-1);
for N = 2:N_Errornodes
[~, u_num, ~, ~, K, ~] = Numsol(E, rho, Ri, Ry, Omega, Poisson, N);
[~, u_exact, ~, ~] = exactsol(E, rho, Ri, Ry, Omega, Poisson, N);
Feeeeel(N-1) = norm(u_num-u_exact')/(N-1);
Cond(N-1) = 1./rcond(K(1:N,1:N));
% http://eprints.ma.man.ac.uk/1997/01/covered/MIMS_ep2013_35.pdf
end
ErrorTable = table((2:N_Errornodes)',Feeeeel', Cond', 'VariableNames',{'Nodes','Error','Condition_number'});
writetable(ErrorTable, 'Errortable.xls')
end