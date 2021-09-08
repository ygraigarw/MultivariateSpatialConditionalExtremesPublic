function Z = RsdSlv(U,D);
% function Z = RsdSlv(U,D);
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Takes uniform random variables (for q variates) and values of delta (for each variate) as input
% Outputs the corresponding STANDARD delta-Laplace (DL(0,1,delta) variates (with mean zero and variance 1)
%
% delta-Laplace parameterised in terms of mean and variance
% Tests show that this function is not accurate for abs(U)>3 !!!!
% You can transform standard DL random variates to DL(mu, sigma^2) simply using "y=x*sigma+mu", 
%  as you would for Gaussians
%
% Input
% U       n x q uniform random values (corresponding to Phi(z_j^N)
% D       q x 1 values of delta for each variate 
%
% Output
% Z       n x q corresponding

options=optimset('Display','off');

[nRls,q]=size(U);

Z=nan(nRls,q);
for iR=1:nRls;
   for j=1:q;
      Kpp=sqrt(gamma(1/D(j))/gamma(3/D(j))); %Kpp scale factor for standard DL 
      if U(iR,j)<0.5; %lower half of distribution
         t = fsolve(@(y)f1(y,U(iR,j),D(j)),1,options);
         Z(iR,j) = -abs(t)^(1/D(j))*Kpp; %standard DL variate DL(0,1,delta)
      elseif U==0.5; %upper half of distribution
         Z(iR,j)=0;
      else;
         t = fsolve(@(y)f2(y,U(iR,j),D(j)),1,options);
         Z(iR,j) = abs(t)^(1/D(j))*Kpp; %standard DL variate DL(0,1,delta)
      end;
      
   end;
end;

function T = f1(y,U,D);
%Note that MATLAB reverses order of arguments in gammainc, and includes gamma scale factor!
%This equation is correct - checked by Phil 20190506
T = 2*U-1+gammainc(y,1/D);

function T = f2(y,U,D);
%Note that MATLAB reverses order of arguments in gammainc, and includes gamma scale factor!
%This equation is correct - checked by Phil 20190506
T = 2*U-1-gammainc(y,1/D);

return;