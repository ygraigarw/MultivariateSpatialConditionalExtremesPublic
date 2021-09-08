function [A,B,M,S,D,R,H]=ABMSDR(X,P,L);
% function [A,B,M,S,D,R,H]=ABMSDR(X,P,L);
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Calculate alpha, beta, mu, sigma, R=(rho1,rho2,lambda) and H for parameters P and distances H

%% HV is a structure containing vector location of reference location (G.Rfr)
%% and remote locations (G.Rmt)
% X.Rfr 1 x 2
% X.Rmt q x 2 (q is number of remote locations)
q=size(X.Rmt,1);
H.H0=nan(q,1);

for j=1:q;
    %% Use local Cartesian coordinates
    tAvrLtt=(X.Rmt(j,2)+X.Rfr(2))/2;
    tDst=[cos(tAvrLtt*pi/180)*(X.Rmt(j,1)-X.Rfr(1)),X.Rmt(j,2)-X.Rfr(2)];
    H.H0(j)=sqrt(tDst*tDst');
end;

%% Critical this ordering matches with definition of H 
for j1=1:q;
   for j2=1:q;
       %% Use local Cartesian coordinates
       tAvrLtt=(X.Rmt(j2,2)+X.Rmt(j1,2))/2;
       tDst=[cos(tAvrLtt*pi/180)*(X.Rmt(j2,1)-X.Rmt(j1,1)),X.Rmt(j2,2)-X.Rmt(j1,2)];
       H.HR(j1,j2)=sqrt(tDst*tDst');
   end;
end;

%% Piecewise linear forms for A, B, M and S
tH=H.H0;
tDlt=L.Dlt;
nH=size(tH,1);
A=nan(nH,1);
B=nan(nH,1);
M=nan(nH,1);
S=nan(nH,1);
D=nan(nH,1);
p=q/3; %p is number of remote locations here

for i=1:p;
    tLct=floor(tH(i)/tDlt)+1;
    for k=1:3; %loop over quantities
        if tLct<=L.n; %there are sufficient spline centres
            if rem(tH(i)/tDlt,1)==0;
                A(i+p*(k-1))=P.A(tLct,k);
                B(i+p*(k-1))=P.B(tLct,k);
                M(i+p*(k-1))=P.M(tLct,k);
                S(i+p*(k-1))=P.S(tLct,k);
                D(i+p*(k-1))=P.D(tLct,k);
            elseif tLct+1<=L.n;
                tInc=1-rem(tH(i)/tDlt,1);
                A(i+p*(k-1))=tInc*P.A(tLct,k)+(1-tInc)*P.A(tLct+1,k);
                B(i+p*(k-1))=tInc*P.B(tLct,k)+(1-tInc)*P.B(tLct+1,k);
                M(i+p*(k-1))=tInc*P.M(tLct,k)+(1-tInc)*P.M(tLct+1,k);
                S(i+p*(k-1))=tInc*P.S(tLct,k)+(1-tInc)*P.S(tLct+1,k);
                D(i+p*(k-1))=tInc*P.D(tLct,k)+(1-tInc)*P.D(tLct+1,k);
            else;
                A(i+p*(k-1))=P.A(L.n,k);
                B(i+p*(k-1))=P.B(L.n,k);
                M(i+p*(k-1))=P.M(L.n,k);
                S(i+p*(k-1))=P.S(L.n,k);
                D(i+p*(k-1))=P.D(L.n,k);
            end;
        else; %insufficient spline centres - use constant value
            A(i+p*(k-1))=P.A(L.n,k);
            B(i+p*(k-1))=P.B(L.n,k);
            M(i+p*(k-1))=P.M(L.n,k);
            S(i+p*(k-1))=P.S(L.n,k);
            D(i+p*(k-1))=P.D(L.n,k);
        end;
    end; %over k
end;

R=P.R';

return;