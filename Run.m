function Run(nNep,nNod,iCndLct,nItr,iRpt);
% function Run(nNep,nNod,iCndLct,nItr,iRpt);
%
% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021
%
% Script to run code to generate a sample (or read data file) and make MCMC inference
% - Allows three variables at common set of locations
% - Allows alpha>1 for subasymptotic levels
% - Ignores conditional quantile constraints for Keef et al.
%
%INPUTS
%nNep      %Non-exceedance probability
%nNod      %Number of nodes
%iCndLct   %Conditioning location
%nItr      %Number of iterations
%iRpt      %Reference number for repeat analysis

%% Note: this version does not use a conditioning variable. First variable or column assumed to be conditioning

%% Setup
pRndSed;

%% Transfer input arguments
X.Nep=nNep;
L.n=nNod;
C.nI=nItr;

%% POTENTIAL USER INPUT: Home directory for work
%HomPth='Z:\individual\philip.jonathan\2021_Satellite\GitAnalysis';
HomPth=pwd;

%% Create analysis directory and cd there
AnlTag=sprintf('NEP%gnNod%giCndLct%giRpt%g',100*X.Nep,L.n,iCndLct,iRpt);
AnlPth=fullfile(HomPth,AnlTag);
mkdir(HomPth,AnlTag);
cd(AnlPth);

fprintf(1,'Multivariate spatial conditional extremes %s started at %s\n',AnlTag,datestr(now,30));

X.DatNam=AnlTag;
X.iCndLct=iCndLct;

%% POTENTIAL USER INPUT: Define data to use
% Load data
if 1;
    load(fullfile(HomPth,'GoldMetNorDat'),'MNUnf','Dat');
    fprintf(1,'Analysis conditioned on SatWnd\n');
else;
    load(fullfile(HomPth,'GoldMNCndNorWav'),'MNUnf','Dat');
    fprintf(1,'Analysis conditioned on NorWav\n');
end;
% Always put conditioning location first
X.AllLct=(1:size(MNUnf,3))';
IsCnd=ismember(X.AllLct,iCndLct);
PltInd=[iCndLct;X.AllLct(IsCnd==0)];
% Manage which variable to condition on, and put variables in right order
Unf1=permute(MNUnf(:,1,PltInd),[1 3 2]);
Unf2=permute(MNUnf(:,2,PltInd),[1 3 2]);
Unf3=permute(MNUnf(:,3,PltInd),[1 3 2]);
%Put LngLtt in the right order
Lct=Dat.LngLtt(PltInd,:);
% Create data set for SCE analysis
DatU=[Unf1 Unf2 Unf3]; % n x 3p array
X.L0=[Lct;Lct;Lct];    % 3p locations on the x-axis in 2D space
X.DatL=-sign(DatU-0.5).*log(1 - 2.*abs(DatU-0.5)); %Uniform to Laplace
X.p=size(X.DatL,2)/3;  %number of different locations
X.q=size(X.DatL,2)-3;  %number of conditioned variables (ie 3(p-1))
% Co-ordinates of reference and remote locations
X.Rfr=X.L0(1,:); %reference location
X.Rmt=X.L0([2:X.p X.p+2:2*X.p 2*X.p+2:3*X.p],:); %3(p-1) "remote" locations, including variables at the conditioning location
X.iVrb=[ones(X.p-1,1);2*ones(X.p-1,1);3*ones(X.p-1,1)]; %variable identifier (conditioning is always done on variable 1)
% Final data for analysis
X.ThrLpl=-sign(X.Nep-0.5).*log(1 - 2.*abs(X.Nep-0.5));
tKep=X.DatL(:,1)>X.ThrLpl;
X.X0C=X.DatL(tKep,1);
X.XRC=X.DatL(tKep,[2:X.p X.p+2:2*X.p 2*X.p+2:3*X.p]);

%% Plot locations of data
set(groot,'DefaultAxesColorOrder',hsv(X.q));
figure(1); clf; hold on;
plot(X.L0(:,1),X.L0(:,2),'r.','markersize',40);
for j=1:size(Lct,1);
    text(Lct(j,1),Lct(j,2),sprintf('%g',j));
end;
pGI(sprintf('%s-Locations%g%s',X.DatNam,100*X.Nep,datestr(now,30)),2);

%% MCMC Setup
C.NgtStr=0.1; %candidate random walk standard deviation
C.NmbToPlt=1000; %number of iterations from end of chain to plot
C.AdpItr1=250; % 100 is the minimum value for this
C.AdpBet=0.05; %Roberts-Rosenthaal adaptive beta

%% Limits for parameters
Lmt.A=[-0.1 1.1]; %limit for each quantity and distance
Lmt.B=[-0.1 1]; 
Lmt.M=[-1 1];
Lmt.S=[0 sqrt(2)+0.5];
Lmt.D=[0.5 3.0]; 
%
Lmt.R=[0 1]; %rho and kappa
Lmt.L=[0 1]; %lambda

%% Handle variable number of M and D parameters
L.HMxm=15; %this is in units of EarthRadius*pi/180 (so approximately degrees)
L.Dlt=L.HMxm/(L.n-1); %distance increment

%% Specify starting solution
% Starting parameters
% alpha, beta, mu, sigma, delta for quantity 1
% alpha, beta, mu, sigma, delta for quantity 2
% alpha, beta, mu, sigma, delta for quantity 3
% 15 rho, kappa, lambda parameters (see ordering list in SceCrr)
L.Prm0=[...
    ones(L.n,1)*0.7;...
    ones(L.n,1)*0.2;...
    ones(L.n,1)*0.0;...
    ones(L.n,1)*0.5;...
    ones(L.n,1)*1.1;...
    ones(L.n,1)*0.7;...
    ones(L.n,1)*0.2;...
    ones(L.n,1)*0.0;...
    ones(L.n,1)*0.5;...
    ones(L.n,1)*1.1;...
    ones(L.n,1)*0.7;...
    ones(L.n,1)*0.2;...
    ones(L.n,1)*0.0;...
    ones(L.n,1)*0.5;...
    ones(L.n,1)*1.1;...
    0.1;...
    0.7*2/5;...
    0.1;...
    0.9*2/5;...
    0.1;...
    0.9*2/5;...
    0.5;...
    0.5*2/5;...
    0.5;...
    0.5*2/5;...
    0.5;...
    0.5*2/5;...
    0.0;...
    0.0;...
    0.0;...
    ];

%% Run MCMC code
MCMC(X,Lmt,C,L);

cd(HomPth);

return;