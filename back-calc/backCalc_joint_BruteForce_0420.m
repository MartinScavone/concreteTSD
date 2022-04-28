function [c,LTE,finalSSE,domainLTE,domainC,domainSSE] = backCalc_joint_BruteForce_0420(TSD_points,TSDvy,TSDvx,estimateC,h,k,G,E,nu,LOAD,pressure,verbose_solve)
% function [c,LTE,finalSSE,domainLTE,domainC,domainSSE] =  backCalc_joint_BruteForce_0420(TSD_points,TSDvy,TSDvx,estimateC,h,k,G,E,nu,LOAD,pressure,verbose_solve)
%
%Back-calcualtion of the joint's LTE index and exact location from nearby
%TSD deflection vertical velocity measurements%
%
%INPUT:
%	TSD_points - location of the TSD sensors where the defl. slope values are taken [m]
%   TSDvy - vertical deflection velocity measurements at the joint's vicinity [m/sec]
%   TSDvx - TSD traveling speed at the time of measurement [m/sec]
%   h - concrete slab thickness [m]
%   k - subgrade's modulus of reaction [N/m3 = 1e9* MPa/mm]
%   G - subgrade's shear strength [N/m2] - for Winkler foundation, G = 0
%   E = concrete slab elastic modulus [N/m2]
%   nu= concrete's Poisson coef. Dimless, default 0.19-0.23
%   LOAD=applied load by the TSD 1/2 axle [N]
%   Pressure-TSD tire pressure [N/m2]
%   verboseness: 1 if you want the solver to be verbose, 0 otherwise
%
%OUTPUT
%   c = exact location of the joint from the TSD axle [m]
%   LTE= joint's LTE index.
%   SSE= [optional] Error domain SSE(domainC,domainLTE)
%
%%Use Van Cauwelaert's 2004 formulation for the jointed slab deflection problem
%
%release candidate v2022-05-01
 
%% preparation 
%define the TSD points to get the TSD deflection slopes.
TSD_points = TSD_points(:);
TSDvy = TSDvy(:);
yDomain = 0;

%% Brute-force back-calculation. Use the [domainC,domainLTE] to compute SSE 
%and the gradient for each back-calculation (because it's so bugging that
%it never converges and climbs the gradient...

%update v2022-0309 -> use the shrank domainC as given by estimateC
domainC = estimateC;
domainLTE = 0.06:0.01:1;
domainSSE = zeros(length(domainC),length(domainLTE));
if verbose_solve
    disp('computing domain for SSE')
    fprintf('\t joint location is bounded between %g and %g \n',max(estimateC),min(estimateC))
end

%% launch brute-force search
% Implementation v2022-04-21
for i = 1:length(domainC)
    if verbose_solve 
        fprintf('\t progress %g percent \n',100.*i/length(domainC))
    end
    for j = 1:length(domainLTE)        
        % simulate vy for this pair of c, LTE
        simulatedVY = getDeflectionVelocity_singleProfile(TSD_points,yDomain,TSDvx,h,E,nu,k,G,LOAD,pressure,domainC(i),domainLTE(j));
        simulatedVY = simulatedVY(:);     %output in m/sec
        
        %SSE term [compute for TSDvy and simulated VY in mm/sec]
        auxSSE = (1e3.*TSDvy - 1e3.*simulatedVY).^2;
        domainSSE(i,j) = sum(auxSSE);
    end
end

%% Locate the minimum of SSE and report its location
%Source -- https://www.mathworks.com/matlabcentral/answers/203754-minimum-value-row-and-column
[argminByRow,argminByCol] = find(domainSSE == min(min(domainSSE)));
c = domainC(argminByRow);
LTE = domainLTE(argminByCol);
finalSSE = min(min(domainSSE));


