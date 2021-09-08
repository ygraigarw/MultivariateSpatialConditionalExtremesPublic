function Nll = SceNll(P,L,X)
% function Nll = SceNll(P,L,X)
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Evaluate negative log-likelihood of MSCE model

%% Evaluate parameters at obervations
[A,B,M,S,D,~,H]=ABMSDR(X,P,L); %all parameters are q x 1

%% Set-up
[n,q]=size(X.XRC); %number of observations and remote locations
PdfZ=nan(n,q);CdfZ=nan(n,q);QntNrm=nan(n,q);PdfNrm=nan(n,q); %initialise arrays

%% Evaluate correlation matrix
Crr=SceCrr(X,P,H);

invCrr=inv(Crr);
invCrrSR=sqrtm(invCrr);
lgrDtrCrr=log(abs(det(Crr)));
if isreal(invCrrSR)==0;
    Nll=inf;
    return;
end;

%% Pdf and cdf
for j=1:q;
    Kpp=sqrt(gamma(1/D(j))/gamma(3/D(j)));
    tM=X.X0C*A'+(X.X0C*ones(1,q)).^(ones(n,1)*B').*(ones(n,1)*M');
    tS=(X.X0C*ones(1,q)).^(ones(n,1)*B').*(ones(n,1)*S')*Kpp;
    PdfZ(:,j)=(D(j)./(2*tS(:,j).*gamma(1/D(j)))).*exp(-abs((X.XRC(:,j)-tM(:,j))./tS(:,j)).^D(j));
    CdfZ(:,j)=0.5+0.5*sign(X.XRC(:,j)-tM(:,j)).*gammainc(abs((X.XRC(:,j)-tM(:,j))./tS(:,j)).^D(j),1/D(j));
    QntNrm(:,j)=norminv(CdfZ(:,j));
    PdfNrm(:,j)=normpdf(QntNrm(:,j));
end; %gammainc as defined by MATLAB has /Gamma; be careful - see doc

%% Evaluate the negative log likelihood
tnll=nan(n,1);
for i=1:n;
    tnll(i)=(q/2)*log(2*pi)+(1/2)*lgrDtrCrr+(1/2)*sum((invCrrSR*QntNrm(i,:)').^2)-sum(log(PdfZ(i,:)./PdfNrm(i,:))); %need to check this is faster
end;
Nll=sum(tnll);
if isnan(Nll)==1;
    Nll=inf;
end;

return;