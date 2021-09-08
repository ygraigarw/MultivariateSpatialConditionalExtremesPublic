function IsOK=PrmChc(Prm,PrmLmt);
% function IsOK=PrmChc(Prm,PrmLmt);
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Check that parameter values are valid (ie within limits)
% All input as arrays not structures

IsOK=0;

t1=Prm<PrmLmt(:,1);
t2=Prm>PrmLmt(:,2);

if sum(t1+t2)==0;
   IsOK=1;
end;

return;