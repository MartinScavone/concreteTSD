function [gradJ,HessianJ] = computeGradient_ContinuousVy(xDomain,yDomain,TSD_meas,h,nu,E,G,k,LOAD,pressure,vx)
%function [gradJ,HessianJ] = computeGradient_ContinuousVy(xDomain,yDomain,TSD_meas,h,nu,E,G,k,LOAD,pressure,vx)
%
%Auxiliary function to the back-calculation tool. 
%Compute the gradient of the cost function for the deflection velocity gradient descent.
%THIS FUNCTION COMPUTES THE GRADIENT FOR THE CONTINUOUS SLAB ONLY (AHEAD OF
%THE JOINT, WHERE I'M OPTIMIZING FOR E,k,G).
%The cost function to minimize is SSE - [TSD_meas - vy(k,E,G)]^2
%
%PROCEED NUMERICALLY, using centered finite differences: 
%Basically, Compute the cost function for the given case +/-1% increase 
%to each variable one at a time to get the partial derivatives over 
%the variables of interest.
%
%INPUT
%   xDomain:    caculation domain, longitudinal direction [m]
%   yDomain:    calculation domain, transverse direction [m]
%   TSD_meas:   Measured deflection velocity at (xDomain, yDOmain) [m/sec]
%   h:          slab thickness [m]
%   E:          slab elastic modulus [N/m2]
%   nu:         slab Poisson coefficient [dimless]
%   k:          subgrade's modulus of reaction [N/m3]
%   G:          subgrade's shear modulus [N/m2]. For Winkler foundation, G =0
%   LOAD:       amount of load [newtons]
%   pressure:   load pressure [N/m2]
%   vx:         TSD travel speed [scalar, m/sec]  
%
%OUTPUT:
%   gradJ = (approximated) gradient of the cost function for the 6
%   variables, at their current values      
%      gradJ[1]: partialJ/partial_k [subgrade's modulus]
%      gradJ[2]: partialJ/partial_E [Concrete Modulus]
%      gradJ[3]: partialJ/partial_g [dim-less quantity related to Pasternak's G - Van C. chap 15]
%
%       HessianJ[1] = partial2J/partialk2
%       HessianJ[2] = partial2J/partialE2
%       HessianJ[3] = partial2J/partialg2
%
%release candidate V2022-05-01

%% CODE BEGINS
% pass TSD_meas to mm/sec (because my SSE is all throughout based on vy in
% mm/sec]
TSD_meas = TSD_meas.*1e3;

%Initialize output
gradJ = zeros(3,1);
HessianJ = zeros(3,1);

%% compute base case [needed for HJJ]

%1) compute vy_base
vy_base = getDeflectionVelocity_continuousSlab(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure);
vy_base = vy_base.*1e3;  %pass from m/sec to mm/sec

%2) compute SSE base
SSE_base = vy_base - TSD_meas;
SSE_base = sum(SSE_base.^2);

%% compute the cost for the ALTERED k scenario
deltaK = 0.01.*k;

%1) Compute vy_plus
vy_plus = getDeflectionVelocity_continuousSlab(xDomain,yDomain,vx,h,E,nu,k+deltaK,G,LOAD,pressure);
vy_plus = vy_plus.*1e3;  %pass from m/sec to mm/sec
%2) compute SSE_plus
SSE_plus = vy_plus - TSD_meas;
SSE_plus = sum(SSE_plus.^2);

%3) compute vy_minus
vy_minus = getDeflectionVelocity_continuousSlab(xDomain,yDomain,vx,h,E,nu,k-deltaK,G,LOAD,pressure);
vy_minus = vy_minus.*1e3;  %pass from m/sec to mm/sec
%4) compute SSE_minus
SSE_minus = vy_minus - TSD_meas;
SSE_minus = sum(SSE_minus.^2);

%5) compute gradJ(1) and HessianJ(1)
gradJ(1) = (SSE_plus - SSE_minus)./(2.*deltaK);
HessianJ(1) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaK.^2);

%% compute the cost for the ALTERED E scenario
deltaE = 0.01.*E;

%1) Compute vy_plus
vy_plus = getDeflectionVelocity_continuousSlab(xDomain,yDomain,vx,h,E+deltaE,nu,k,G,LOAD,pressure);
vy_plus = vy_plus.*1e3;   %pass from m/sec to mm/sec
%2) compute SSE_plus
SSE_plus = vy_plus - TSD_meas;
SSE_plus = sum(SSE_plus.^2);

%3) compute vy_minus
vy_minus = getDeflectionVelocity_continuousSlab(xDomain,yDomain,vx,h,E-deltaE,nu,k,G,LOAD,pressure);
vy_minus = vy_minus.*1e3;   %pass from m/sec to mm/sec
%4) compute SSE_minus
SSE_minus = vy_minus - TSD_meas;
SSE_minus = sum(SSE_minus.^2);

%5) compute gradJ(1) and HessianJ(1)
gradJ(2) = (SSE_plus - SSE_minus)./(2.*deltaE);
HessianJ(2) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaE.^2);

%% compute the cost for the ALTERED G scenario 
%BUG FIX 2021-06-15: If WINKLER FOUNDATION (G = 0), bypass this operation
%and return a zero, cos IT WOULD OTHERWISE MAKE DIVISION BY ZERO!

if G>0   
    deltaG = 0.01.*G;

    %1) Compute vy_plus
    vy_plus = getDeflectionVelocity_continuousSlab(xDomain,yDomain,vx,h,E,nu,k,G+deltaG,LOAD,pressure);
    vy_plus = vy_plus.*1e3;   %pass from m/sec to mm/sec
    %2) compute SSE_plus
    SSE_plus = vy_plus - TSD_meas;
    SSE_plus = sum(SSE_plus.^2);

    %3) compute vy_minus
    vy_minus = getDeflectionVelocity_continuousSlab(xDomain,yDomain,vx,h,E,nu,k,G-deltaG,LOAD,pressure);
    vy_minus = vy_minus.*1e3;   %pass from m/sec to mm/sec
    %4) compute SSE_minus
    SSE_minus = vy_minus - TSD_meas;
    SSE_minus = sum(SSE_minus.^2);

    %5) compute gradJ(1) and HessianJ(1)
    gradJ(3) = (SSE_plus - SSE_minus)./(2.*deltaG);
    HessianJ(3) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaG.^2);
    

else
    gradJ(3) = 0;
    HessianJ(3) = 0;
end        
end %<--- endfunction