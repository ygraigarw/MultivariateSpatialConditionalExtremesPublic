function PltDgn(C,iI);
%% function PltDgn(C,iI);
%
% Conditional spatial extremes with delta Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Diagnostic plot of alpha against beta, proposal std and acceptance rates

clf;
q=size(C.Alp,2); %number of remote variables
p=q/3; %number of remote locations
tClr=jet(p);
tMrk={'s';'^';'o'};

%% alpha on beta
subplot(1,2,1);hold on;
jL=0;
for iV=1:3;
    for iL=1:p;
        jL=jL+1;
        plot(C.Alp(max(1,iI-19):iI,jL),C.Bet(max(1,iI-19):iI,jL),'--','color',tClr(iL,:));
        plot(C.Alp(iI,jL),C.Bet(iI,jL),'marker',tMrk{iV},'color',tClr(iL,:),'markersize',10,'markerfacecolor',tClr(iL,:)); title('\beta on \alpha');
    end
end;
pDfl; pAxsLmt;
set(gca,'xlim',[min(C.Alp(iI,:))-0.05,max(C.Alp(iI,:))+0.05]);
set(gca,'ylim',[min(C.Bet(iI,:))-0.05,max(C.Bet(iI,:))+0.05]);
title('\beta on \alpha');

%% proposal std
subplot(2,2,2); 
plot(C.Ngt(1:iI,:)); title('Proposal std'); 
pDfl; pAxsLmt;

%% acceptance rates
subplot(2,2,4); 
plot(C.AccRat(1:iI,:)); title('Acceptance rate'); 
pDfl; pAxsLmt;

drawnow;

return;