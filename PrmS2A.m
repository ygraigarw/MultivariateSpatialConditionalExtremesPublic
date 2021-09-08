function Prm=PrmS2A(L,Q);
% function P=PrmS2A(L,Q)
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Convert alpha, beta, mu, sigma and R=(rho1,rho2,lambda) from structure to array

n=L.n;

%% Ordering of parameters used is
% A, B, M and S for quantity 1, then
% A, B, M and S for quantity 2, then
% A, B, M and S for quantity 3, then
% R (ie the common parameters rho1, rho2, lambda)

if size(Q.A,1)>1; %must be a column vector, so P
   Prm=[Q.A(:,1);Q.B(:,1);Q.M(:,1);Q.S(:,1);Q.D(:,1);Q.A(:,2);Q.B(:,2);Q.M(:,2);Q.S(:,2);Q.D(:,2);Q.A(:,3);Q.B(:,3);Q.M(:,3);Q.S(:,3);Q.D(:,3);Q.R'];
else; %must be Lmt
   Prm=[ones(n,1)*Q.A;ones(n,1)*Q.B;ones(n,1)*Q.M;ones(n,1)*Q.S;ones(n,1)*Q.D;ones(n,1)*Q.A;ones(n,1)*Q.B;ones(n,1)*Q.M;ones(n,1)*Q.S;ones(n,1)*Q.D;ones(n,1)*Q.A;ones(n,1)*Q.B;ones(n,1)*Q.M;ones(n,1)*Q.S;ones(n,1)*Q.D;ones(15,1)*Q.R];
end;
   
return;