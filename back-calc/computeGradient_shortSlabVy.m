function [gradJ,HessianJ] = computeGradient_shortSlabVy(xDomain,yDomain,TSD_meas,h,nu,E,G,k,c,LTE,LossSupport,LOAD,pressure,vx)
%function [gradJ,HessianJ] = computeGradient_shortSlabVy(xDomain,yDomain,TSD_meas,h,nu,E,G,k,c,LTE,LossSupport,LOAD,pressure,vx)
%
%Auxiliary function to the back-calculation tool. 
%Compute the gradient of the cost function for the deflection velocity gradient descent.
%THIS FUNCTION COMPUTES THE GRADIENT FOR THE SHORT SLAB CASE, WHERE I'M
%BACK-CALCUALTING SIMULTANEOUSLY FOR k,E,G, and the joint's LTE.
%The cost function to minimize is SSE - [TSD_meas - vy(k,E,G,LTE|c,h,nu)]^2
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
%   G:          subgrade's shear modulus [N/m]. For Winkler foundation, G =0
%   c:          distance between load center and joint [m]
%   LTE:        load transfer efficiency index at the joint (deflection-based) [0-1 variable].
%   LossSupport: Subgrade's loss-of-support index [0-1 variable]
%   LOAD:       amount of load [newtons]
%   pressure:   load pressure [N/m2]
%   vx:         TSD travel speed [scalar, m/sec]  
%
%OUTPUT:
%   gradJ = (approximated) gradient of the cost function for the target
%   variables, at their current values   
%      gradJ[1]: partialJ/partial_k [subgrade's modulus]
%      gradJ[2]: partialJ/partial_E [Concrete Modulus]
%      gradJ[3]: partialJ/partial_G [subgrade's linear shear strength modulus]
%      gradJ[4]: partialJ/partial_LTE

%       HessianJ[1] = partial2J/partialk2
%       HessianJ[2] = partial2J/partialE2
%       HessianJ[3] = partial2J/partialG2
%       HessianJ[4] = partial2J/partial_LTE2
%
%V1.0 - 2023-04-25 - Giulietta or not
%   First version of the code, based on back-calculation solver 2023-01-25

%% CODE BEGINS
% pass TSD_meas to mm/sec (because my SSE is all throughout based on vy in
% mm/sec]
TSD_meas = TSD_meas.*1e3;

%Initialize output
gradJ = zeros(4,1);
HessianJ = zeros(4,1);

%% compute base case [needed for HJJ]

%1) compute vy_base
vy_base = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LossSupport,LTE);
vy_base = vy_base.*1e3;  %pass from m/sec to mm/sec

%2) compute SSE base
SSE_base = vy_base - TSD_meas;
SSE_base = sum(SSE_base.^2);


%% compute the cost for the ALTERED k scenario
deltaK = 0.01.*k;

%1) Compute vy_plus
vy_plus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k+deltaK,G,LOAD,pressure,c,LossSupport,LTE);
vy_plus = vy_plus.*1e3;   %pass from m/sec to mm/sec
%2) compute SSE_plus
SSE_plus = vy_plus - TSD_meas;
SSE_plus = sum(SSE_plus.^2);

%3) compute vy_minus
vy_minus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k-deltaK,G,LOAD,pressure,c,LossSupport,LTE);
vy_minus = vy_minus.*1e3;   %pass from m/sec to mm/sec
%4) compute SSE_minus
SSE_minus = vy_minus - TSD_meas;
SSE_minus = sum(SSE_minus.^2);

%5) compute gradJ(1) and HessianJ(1)
gradJ(1) = (SSE_plus - SSE_minus)./(2.*deltaK);
HessianJ(1) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaK.^2);


%% compute the cost for the ALTERED E scenario
deltaE = 0.01.*E;

%1) Compute vy_plus
vy_plus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E+deltaE,nu,k,G,LOAD,pressure,c,LossSupport,LTE);
vy_plus = vy_plus.*1e3;   %pass from m/sec to mm/sec
%2) compute SSE_plus
SSE_plus = vy_plus - TSD_meas;
SSE_plus = sum(SSE_plus.^2);

%3) compute vy_minus
vy_minus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E-deltaE,nu,k,G,LOAD,pressure,c,LossSupport,LTE);
vy_minus = vy_minus.*1e3;   %pass from m/sec to mm/sec
%4) compute SSE_minus
SSE_minus = vy_minus - TSD_meas;
SSE_minus = sum(SSE_minus.^2);

%5) compute gradJ(2) and HessianJ(2)
gradJ(2) = (SSE_plus - SSE_minus)./(2.*deltaE);
HessianJ(2) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaE.^2);


%% compute the cost for the ALTERED G scenario
%BUG FIX 2021-06-15: If WINKLER FOUNDATION (G = 0), bypass this operation
%and return a zero, cos IT WOULD OTHERWISE MAKE DIVISION BY ZERO!

if G>0   
    deltaG = 0.01.*G;
    %1) Compute vy_plus
    vy_plus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k,G+deltaG,LOAD,pressure,c,LossSupport,LTE);
    vy_plus = vy_plus.*1e3;   %pass from m/sec to mm/sec
    %2) compute SSE_plus
    SSE_plus = vy_plus - TSD_meas;
    SSE_plus = sum(SSE_plus.^2);

    %3) compute vy_minus
    vy_minus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k,G-deltaG,LOAD,pressure,c,LossSupport,LTE);
    vy_minus = vy_minus.*1e3;   %pass from m/sec to mm/sec
    %4) compute SSE_minus
    SSE_minus = vy_minus - TSD_meas;
    SSE_minus = sum(SSE_minus.^2);

    %5) compute gradJ(3) and HessianJ(3)
    gradJ(3) = (SSE_plus - SSE_minus)./(2.*deltaG);
    HessianJ(3) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaG.^2);
else
    gradJ(3) = 0;
    HessianJ(3) = 0;
end

%% compute the cost for the ALTERED LTE scenario
deltaLTE = 0.01.*LTE;

%1) Compute vy_plus
vy_plus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LossSupport,LTE+deltaLTE);
vy_plus = vy_plus.*1e3;   %pass from m/sec to mm/sec
%2) compute SSE_plus
SSE_plus = vy_plus - TSD_meas;
SSE_plus = sum(SSE_plus.^2);

%3) compute vy_minus
vy_minus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LossSupport,LTE-deltaLTE);
vy_minus = vy_minus.*1e3;   %pass from m/sec to mm/sec
%4) compute SSE_minus
SSE_minus = vy_minus - TSD_meas;
SSE_minus = sum(SSE_minus.^2);

%5) compute gradJ(4) and HessianJ(4)
gradJ(4) = (SSE_plus - SSE_minus)./(2.*deltaLTE);
HessianJ(4) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaLTE.^2);


end %<--- endfunction