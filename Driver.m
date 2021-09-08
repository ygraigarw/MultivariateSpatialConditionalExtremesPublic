% MULTIVARIATE SPATIAL CONDITIONAL EXTREMES
% Philip Jonathan, Rob Shooter, Emma Ross
% September 2021

%% Plotting utilities available from GitHub at https://github.com/ygraigarw/pMtlUtl
%addpath Y:\individual\Philip.Jonathan\PhilipGit\Utilities;
addpath 'C:\PhilipGit\pMtlUtl';

%% Notes before running
%PLOTTING UTILITIES:
%- These need to be available (above) for the code to run
%HOME DIRECTORY:
%- You can optionally change home directory manually in Run.m
%- As specified, all code and data is expected in current directory
%DATA:
%- Code expects three quantities for MSCE modelling
%- Sample must be on standard Laplace margins
%- Example data set used in OE paper provided with code (2MB)
%SCALING OF RESIDUAL CORRELATION PARAMETERS:
%- In SceCrr.m, the values of Scl1 and Scl2 can be adjusted to ensure that values of rho, kappa lambda are in [0,1] for neatness
%SPATIAL SCALE
%- You might want to change L.HMxm in Run.m. This controls the total distance over which distance parameters are estimated


%% Specify run parameters
Nep=0.95;  % Non-exceedance probability tau
nNod=6;    % n_Nds
iCndLct=1; % conditioning location (1-14)
nItr=1000; % number of MCMC iterations (needs to be >=10k for serious run)
iRpt=1001; % counter number for analysis (any positive integer)

%% Run MSCE
Run(Nep,nNod,iCndLct,nItr,iRpt);