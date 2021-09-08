function [P,PNms]=PrmA2S(L,Prm);
% function P=DAPrmA2S(Prm)
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Convert alpha, beta, mu, sigma, R=(rho1,rho2,lambda) from array to structure

n=L.n;

PNms=cell(10*n+3,1);

j=0;
% Parameters for first quantity
for i1=1:n;
    j=j+1;
    P.A(i1,1)=Prm(j);
    PNms{j}=sprintf('A1(%g)',i1);
end;
for i1=1:n;
    j=j+1;
    P.B(i1,1)=Prm(j);
    PNms{j}=sprintf('B1(%g)',i1);
end;
for i1=1:n
    j=j+1;
    P.M(i1,1)=Prm(j);
    PNms{j}=sprintf('M1(%g)',i1);
end;
for i1=1:n
    j=j+1;
    P.S(i1,1)=Prm(j);
    PNms{j}=sprintf('S1(%g)',i1);
end;
for i1=1:n
    j=j+1;
    P.D(i1,1)=Prm(j);
    PNms{j}=sprintf('D1(%g)',i1);
end;
% Parameters for second quantity
for i1=1:n;
    j=j+1;
    P.A(i1,2)=Prm(j);
    PNms{j}=sprintf('A2(%g)',i1);
end;
for i1=1:n;
    j=j+1;
    P.B(i1,2)=Prm(j);
    PNms{j}=sprintf('B2(%g)',i1);
end;
for i1=1:n
    j=j+1;
    P.M(i1,2)=Prm(j);
    PNms{j}=sprintf('M2(%g)',i1);
end;
for i1=1:n
    j=j+1;
    P.S(i1,2)=Prm(j);
    PNms{j}=sprintf('S2(%g)',i1);
end;
for i1=1:n
    j=j+1;
    P.D(i1,2)=Prm(j);
    PNms{j}=sprintf('D2(%g)',i1);
end;
% Parameters for third quantity
for i1=1:n;
    j=j+1;
    P.A(i1,3)=Prm(j);
    PNms{j}=sprintf('A3(%g)',i1);
end;
for i1=1:n;
    j=j+1;
    P.B(i1,3)=Prm(j);
    PNms{j}=sprintf('B3(%g)',i1);
end;
for i1=1:n
    j=j+1;
    P.M(i1,3)=Prm(j);
    PNms{j}=sprintf('M3(%g)',i1);
end;
for i1=1:n
    j=j+1;
    P.S(i1,3)=Prm(j);
    PNms{j}=sprintf('S3(%g)',i1);
end;
for i1=1:n
    j=j+1;
    P.D(i1,3)=Prm(j);
    PNms{j}=sprintf('D3(%g)',i1);
end;
% Common parameters
for i1=1:15
    j=j+1;
    P.R(i1)=Prm(j);
    PNms{j}=sprintf('R(%g)',i1);
end;

return;