function C=MCMC(X,Lmt,C,L);
% function C=MCMC(X,Lmt,C,L);
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Main MCMC function

figure(2); clf;

%% Specify PrmLmt array
q=size(X.Rmt,1); %number of conditioned variables
PrmLmt=PrmS2A(L,Lmt); %parameter limits as array
nP=size(PrmLmt,1); %number of parameters

%% Local copies of parameters for neatness
nI=C.nI;
AdpItr1=C.AdpItr1; %number of iterations for fixed nuggest
AdpBet=C.AdpBet;
NmbToPlt=C.NmbToPlt;
NgtStr=C.NgtStr;

fprintf(1,'Starting MCMC analysis for case NEP = %g.\n',X.Nep);

%% Allow continuation of completed previous analysis
tFil=sprintf('%s-Prm%.2f.mat',X.DatNam,100*X.Nep);
if exist(fullfile(cd,tFil),'file')==2
    fprintf(1,'Using previous chain\n');
    load(tFil,'C');
    tC=C; clear C;
    
    %Make sure mI is correct
    mI=size(tC.Alp,1);
    tC.AccRat=tC.AccRat(1:mI,:);
    tC.Ngt=tC.Ngt(1:mI,:);
    tC.Nll=tC.Nll(1:mI,:);
    tC.Prm=tC.Prm(1:mI,:);
    
    %MCMC arrays
    C.AccRat=nan(mI+nI,nP);C.AccRat(1:mI,:)=tC.AccRat;
    C.Ngt=nan(mI+nI,nP);C.Ngt(1:mI,:)=tC.Ngt;
    C.Nll=nan(mI+nI,1);C.Nll(1:mI,:)=tC.Nll;
    C.Prm=nan(mI+nI,nP);C.Prm(1:mI,:)=tC.Prm;
    C.nI=mI+nI;
    Nll=C.Nll(mI);
    Prm=C.Prm(mI,:)';
else
    mI=0;
    %MCMC arrays
    C.AccRat=nan(nI,nP);C.AccRat(1,:)=0;
    C.Ngt=nan(nI,nP);
    C.Nll=nan(nI,1);
    C.Prm=nan(nI,nP);
end

%% Loop over MCMC iterations

tic; %count speed

for iI=mI+1:mI+nI;
    
    %Find valid starting solution if iI=1
    if iI==1;
        fprintf(1,'Checking starting solution\n');
        
        %Starting solution is user input
        PrmStr=L.Prm0;
        tP=PrmA2S(L,PrmStr);
        tP.p=q/2+1;
        
        NllStr=SceNll(tP,L,X);
        
        if isinf(NllStr)==1; %no valid starting solution found. Terminate.
            fprintf(1,'Warning: invalid starting solution. Terminating.\n');
            return;
        else %make the current state the starting state for MCMC
            Prm=PrmStr;
            Nll=NllStr;
            fprintf(1,'Starting solution found. Starting MCMC\n');
        end
        
    end;

    %% Iteration counter on screen
    if rem(iI,10)==0;
        fprintf(1,'+');
    else
        fprintf(1,'.');
    end;
    if rem(iI,100)==0;
        fprintf(1,'\n');
    end;
    
    %% Loop over parametric forms
    for iP=1:nP; %Metropolis Hastings in Gibbs, one parameter at a time
        
        %% Define candidate in terms of current state
        PrmC=Prm;
        
        n=L.n;
        if iI<=AdpItr1; %fixed nugget
            tNgt=NgtStr;
            PrmC(iP)=PrmC(iP)+randn*tNgt;
        else; %adaptive Metropolis for A B M S
            if iP<=n; %update alpha beta mu and sigma together
                jP=(0:(3*5-1))'*n+iP;
                nJ=size(jP,1);
                t1=real((1-AdpBet)*2.38*sqrtm(cov(C.Prm(max(1,iI-999):iI-1,jP))/nJ)*randn(nJ,1));
                t2=AdpBet*0.1*(randn(nJ,1)/nJ);
                PrmC(jP)=PrmC(jP)+t1+t2;
                tNgt=NaN;
                C.Ngt(iI-1,jP(2:end))=NaN;
                C.AccRat(iI-1,jP(2:end))=NaN;
            elseif iP>(3*5)*n;
                if C.AccRat(iI-1,iP)>0.3;
                    tNgt=C.Ngt(iI-1,iP)*1.05;
                elseif C.AccRat(iI-1,iP)<0.2;
                    tNgt=C.Ngt(iI-1,iP)*0.95;
                else;
                    tNgt=C.Ngt(iI-1,iP);
                end;
                if tNgt>2; %max proposal sd is unity
                    tNgt=2;
                end;
                PrmC(iP)=PrmC(iP)+randn*tNgt;
            end;
        end;
        
        %% Evaluate likelihood at current state and candidate
        PC=PrmA2S(L,PrmC);
                
        IsOK=PrmChc(PrmC,PrmLmt);
        if IsOK==1;
            NllC=SceNll(PC,L,X);
        else;
            NllC=inf;
        end;
        
        if isreal(NllC)==0;
            fprintf(1,'Imaginary NllC\n');
        end
        
        %% MH acceptance step
        if (exp(-NllC+Nll) > rand) && isinf(NllC)==0 && isnan(NllC)==0;
            Prm=PrmC;
            Nll=NllC;
            if iI>1;
                if iI>100; %Only use last 100 iterations to adjust acceptance rate
                    jI=100;
                else;
                    jI=iI;
                end;
                if iI<=AdpItr1; %fixed nugget
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1)+1)/jI;
                else;
                    if iP<=n;
                        C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1)+1)/jI;
                    elseif iP>10*n;
                        C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1)+1)/jI;
                    end;
                end;
            end;
        else;
            if iI>1; %Only use last 100 iterations to adjust acceptance rate
                if iI>100;
                    jI=100;
                else;
                    jI=iI;
                end;
                if iI<=AdpItr1; %fixed nugget
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1))/jI;
                else;
                    if iP<=n;
                        C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1))/jI;
                    elseif iP>10*n;
                        C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1))/jI;
                    end;
                end;
            end;
        end;
        
        %% Save proposal standard deviation
        C.Ngt(iI,iP)=tNgt;
        
    end;
    
    %% Update after complete iteration over variables
    C.Prm(iI,:)=Prm';
    C.Nll(iI)=Nll;
    [tA,tB]=ABMSDR(X,PrmA2S(L,C.Prm(iI,:)'),L);
    C.Alp(iI,:)=tA';
    C.Bet(iI,:)=tB';
    
    if iI<=50;
        if rem(iI,10)==0;
            
            %% Trace plots
            figure(1); clf;
            PltTrc(C,iI,iI,L); %trace plots
            tFil=sprintf('%s-TrcPrg%.2f',X.DatNam,100*X.Nep);
            pGI(tFil,2);
            
        end;
    end;
    
    if rem(iI,100)==0;
        
        %% Trace plots
        figure(1); clf;
        PltTrc(C,iI,iI,L); %trace plots
        tFil=sprintf('%s-TrcPrg%.2f',X.DatNam,100*X.Nep);
        pGI(tFil,2);
        
        %% Constraint diagnostics
        figure(2); clf;
        PltDgn(C,iI);
        
        %% Parameters with distance and direction
        PltPrm(C,iI,NmbToPlt,L);
        
        %% Save end point for start of new run
        PrmPrg=C.Prm(iI,:)';
        NgtPrg=C.Ngt(iI,:)';
        ItrPrg=iI;
        tFil=sprintf('%s-PrmPrg%.2f.mat',X.DatNam,100*X.Nep);
        save(tFil,'PrmPrg','NgtPrg','ItrPrg');
        
        %% Log
        fprintf(1,'Completed iteration %g. %g iterations in %g seconds\n',iI,100,toc);
        tic;
        
    end;
    
    if rem(iI,1000)==0 || iI==mI+nI;
        
        %% Save whole chain
        DatNam=X.DatNam;
        Nep=X.Nep;
        tFil=sprintf('%s-Prm%.2f.mat',DatNam,100*Nep);
        save(tFil,'C','mI','nI','L','NmbToPlt','DatNam','Nep');
        
        figure(1);
        clf; PltTrc(C,iI,iI,L); % save full trace plot
        pGI(sprintf('%s-trace%g%s',DatNam,100*Nep,datestr(now,30)),2);
        
        clf; PltTrc(C,iI,NmbToPlt,L); % save end of trace plot
        pGI(sprintf('%s-trace_end%g%s',DatNam,100*Nep,datestr(now,30)),2);
        
        figure(2);
        clf; PltDgn(C,iI); % diagnostic plot
        pGI(sprintf('%s-Dgn%g%s',DatNam,100*Nep,datestr(now,30)),2);
        
        tStr=cell(3,1);
        tStr{1}=sprintf('%s-Prm%g%s',DatNam,100*Nep,datestr(now,30));
        tStr{2}=sprintf('%s-DrcAlp%g%s',DatNam,100*Nep,datestr(now,30));
        tStr{3}=sprintf('%s-DrcBet%g%s',DatNam,100*Nep,datestr(now,30));
        PltPrm(C,iI,NmbToPlt,L,tStr);
    end;
    
end;

return;