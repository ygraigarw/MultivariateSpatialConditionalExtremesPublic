function PltPrm(C,nI,nToPlt,L,tStr);
%function PltPrm(C,nI,nToPlt,L,tStr);
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Parameter plots (no need to understand further)

%Handle whether function is generating intermediate or final plots
if nargin==4;
    tStr=[];
end;

nDlt=min(nI,nToPlt);

Prm=C.Prm(nI-nDlt+1:nI,:);
medP=PrmA2S(L,median(Prm)');
lowP=PrmA2S(L,quantile(Prm,0.025)');
uppP=PrmA2S(L,quantile(Prm,0.975)');

Dst=(0:L.Dlt:L.HMxm)';

figure(3); clf;

subplot(2,3,1); hold on;
set(gca,'colororder',[0 0 0;1 0.65 0;1 0 0]);
plot(Dst,lowP.A,'--');
plot(Dst,uppP.A,'--');
plot(Dst,medP.A,'-');
pDfl;pAxsLmt;
title '\alpha';

subplot(2,3,2); hold on;
set(gca,'colororder',[0 0 0;1 0.65 0;1 0 0]);
plot(Dst,lowP.B,'--');
plot(Dst,uppP.B,'--');
plot(Dst,medP.B,'-');
pDfl;pAxsLmt;
title '\beta';

subplot(2,3,3); hold on;
set(gca,'colororder',[0 0 0;1 0.65 0;1 0 0]);
plot(Dst,lowP.M,'--');
plot(Dst,uppP.M,'--');
plot(Dst,medP.M,'-');
pDfl;pAxsLmt;
title '\mu';

subplot(2,3,4); hold on;
set(gca,'colororder',[0 0 0;1 0.65 0;1 0 0]);
plot(Dst,lowP.S,'--');
plot(Dst,uppP.S,'--');
plot(Dst,medP.S,'-');
pDfl;pAxsLmt;
title '\sigma';

subplot(2,3,5); hold on;
set(gca,'colororder',[0 0 0;1 0.65 0;1 0 0]);
plot(Dst,lowP.D,'--');
plot(Dst,uppP.D,'--');
plot(Dst,medP.D,'-');
pDfl;pAxsLmt;
title '\delta';

subplot(2,3,6); hold on;
plot(lowP.R,'kv');
plot(uppP.R,'k^');
plot(medP.R,'ko');
pDfl;pAxsLmt;
set(gca,'TickLabelInterpreter','latex');
set(gca,'xtick',[1:15],'xticklabels',{'$\rho_1$','$\rho_2$','$\rho_3$','$\rho_4$','$\rho_5$','$\rho_6$','$\rho_7$','$\rho_8$','$\rho_9$','$\rho_{10}$','$\rho_{11}$','$\rho_{12}$','$\lambda_1$','$\lambda_2$','$\lambda_3$'});
title 'Others';

if isempty(tStr)==0; pGI(tStr{1},2); end;

return;