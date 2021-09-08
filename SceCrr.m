function Crr=SceCrr(X,P,H);
% function Crr=SceCrr(X,P,H);
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Estimate the (conditional Gaussian) residual correlation matrix

q=size(H.H0,1); %number of remote variables

%% Scaling parameters are hard-coded
% These might need to be adjusted by hand for new problem
Scl1=100;
Scl2=5; %This was =2 in original code / analysis

%P.R(1) = scale for Y1|Y1     rho
%P.R(2) = exponent for Y1|Y1  kappa etc. (see paper)
%P.R(3) = scale for Y2|Y2
%P.R(4) = exponent for Y2|Y2
%P.R(5) = scale for Y3|Y3
%P.R(6) = exponent for Y3|Y3
%
%P.R(7) = scale for (Y1,Y2)
%P.R(8) = exponent for (Y1,Y2)
%P.R(9) = scale for (Y1,Y3)
%P.R(10) = exponent for (Y1,Y3)
%P.R(11) = scale for (Y2,Y3)
%P.R(12) = exponent for (Y2,Y3)
%
%P.R(13) = correlation at zero distance for (Y1,Y2) lambda
%P.R(14) = correlation at zero distance for (Y1,Y3) lambda
%P.R(15) = correlation at zero distance for (Y2,Y3) lambda

iVrb=X.iVrb; %note that by definition the conditioning variate is always 1

%% Calculate correlations with reference location
Rho=nan(q,1);
for iC=1:q;
    tH=H.H0(iC);
    switch iVrb(iC);
        case 1;
            Rho(iC)=exp( -( tH / (Scl1*(P.R(1))) )^(Scl2*P.R(2)) );
        case 2;
            Rho(iC)=P.R(13)*exp( -( tH / (Scl1*(P.R(7))) )^(Scl2*P.R(8)) );
        case 3;
            Rho(iC)=P.R(14)*exp( -( tH / (Scl1*(P.R(9))) )^(Scl2*P.R(10)) );
    end;
end;

%% Calculate correlations between remote locations
Crs=nan(q,q);
for iC=1:q;
    for jC=1:q;
        tH=H.HR(iC,jC);
        switch iVrb(iC)
            case 1;
                switch iVrb(jC);
                    case 1;
                        Crs(iC,jC)=exp( -( tH / (Scl1*(P.R(1))) )^(Scl2*P.R(2)) );
                    case 2;
                        Crs(iC,jC)=P.R(13)*exp( -( tH / (Scl1*(P.R(7))) )^(Scl2*P.R(8)) );
                    case 3;
                        Crs(iC,jC)=P.R(14)*exp( -( tH / (Scl1*(P.R(9))) )^(Scl2*P.R(10)) );
                end;
            case 2;
                switch iVrb(jC);
                    case 1;
                        Crs(iC,jC)=P.R(13)*exp( -( tH / (Scl1*(P.R(7))) )^(Scl2*P.R(8)) );
                    case 2;
                        Crs(iC,jC)=exp( -( tH / (Scl1*(P.R(3))) )^(Scl2*P.R(4)) );
                    case 3;
                        Crs(iC,jC)=P.R(15)*exp( -( tH / (Scl1*(P.R(11))) )^(Scl2*P.R(12)) );
                end;
            case 3;
                switch iVrb(jC);
                    case 1;
                        Crs(iC,jC)=P.R(14)*exp( -( tH / (Scl1*(P.R(9))) )^(Scl2*P.R(10)) );
                    case 2;
                        Crs(iC,jC)=P.R(15)*exp( -( tH / (Scl1*(P.R(11))) )^(Scl2*P.R(12)) );
                    case 3;
                        Crs(iC,jC)=exp( -( tH / (Scl1*(P.R(5))) )^(Scl2*P.R(6)) );
                end;
        end;
    end;
end;
%
Crr=nan(q,q);
for iC=1:q;
    for jC=1:q;
        t1=sqrt(1-Rho(iC)^2);
        t2=sqrt(1-Rho(jC)^2);
        if abs(t1*t2)>1e-10;
            Crr(iC,jC)=(Crs(iC,jC)-Rho(iC)*Rho(jC))/(t1*t2);
        else;
            Crr(iC,jC)=0;
        end;
    end;
end;

return;