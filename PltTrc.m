function PltTrc(C,iI,nToPlt,L);
%% function PltTrc(C,iI,nToPlt);
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% MCMC trace plot (no need to understand further)

nPlt=size(C.Prm,2)+3;
[~,PNms]=PrmA2S(L,C.Prm(iI,:)');

strI=max(0,iI-nToPlt)+1;

if rem(sqrt(nPlt),1)==0
   td1=sqrt(nPlt);
   td2=sqrt(nPlt);
else
   td1=floor(sqrt(nPlt));
   td2=floor(sqrt(nPlt))+2;
end;
for j=1:size(C.Prm,2);
   subplot(td1,td2,j); hold on;
   plot(C.Prm(strI:iI,j)); title(PNms{j}); pDfl; pAxsLmt;
end;
subplot(td1,td2,nPlt-2); plot((1:iI)',C.Nll(1:iI)); title('NLL'); pDfl; pAxsLmt;
subplot(td1,td2,nPlt-1); plot((strI:iI)',C.Nll(strI:iI)); title('NLL'); pDfl; pAxsLmt;
subplot(td1,td2,nPlt); plot((iI-9:iI)',C.Nll(iI-9:iI)); title('NLL'); pDfl; pAxsLmt;
drawnow;

return;