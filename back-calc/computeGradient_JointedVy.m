function [gradJ,HessianJ] = computeGradient_JointedVy(xDomain,yDomain,TSD_meas,h,nu,E,G,k,c,LTE,LossSupport,LOAD,pressure,vx,doLSPT)
%function [gradJ,HessianJ] = computeGradient_JointedVy(xDomain,yDomain,TSD_meas,h,nu,E,G,k,c,LTE,LossSupport,LOAD,pressure,vx,doLSPT)
%
%Auxiliary function to the back-calculation tool. 
%Compute the gradient of the cost function for the deflection velocity gradient descent.
%THIS FUNCTION COMPUTES THE GRADIENT FOR THE JOINTED SLAB, WHERE I'M
%BACK-CALCUALTING FOR Loss of Support & LTE.
%The cost function to minimize is SSE - [TSD_meas - vy(Loss,LTE|c,k,E,G)]^2
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
%   doLSPT:     Boolean: 1 to solver for loss of support at the subgrade, 0 if not.
%
%OUTPUT:
%   gradJ = (approximated) gradient of the cost function for the target
%   variables, at their current values     
%      gradJ[1]: partialJ/partial_LTE
%      gradJ[2]: partialJ/partial_LossSupport
%     
%       HessianJ[1] = partial2J/partial_LTE2
%       HessianJ[2] = partial2J/partial_LossSupport2
%

%V2.1 - 2023-0405 Spring
% Add a boolean input variable "doLSPT" to control whether or not to
% compute Loss of support (and roll back to original back-calc formulation,
% w/o LSPT)

%V2.0 - 2023-01-25 Flu.
% Experimental version to now solve for the loss of support of the
% subgrade and the joint LTE.

%V1.0 Maconha Alleman Aftermath 2022-04-21.
%First build of back-calc function based on vy.

%% CODE BEGINS
% pass TSD_meas to mm/sec (because my SSE is all throughout based on vy in
% mm/sec]
TSD_meas = TSD_meas.*1e3;

%Initialize output
gradJ = zeros(2,1);
HessianJ = zeros(2,1);

%% compute base case [needed for HJJ]

%1) compute vy_base
vy_base = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LossSupport,LTE);
vy_base = vy_base.*1e3;  %pass from m/sec to mm/sec

%2) compute SSE base
SSE_base = vy_base - TSD_meas;
SSE_base = sum(SSE_base.^2);

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

%5) compute gradJ(1) and HessianJ(1)
gradJ(1) = (SSE_plus - SSE_minus)./(2.*deltaLTE);
HessianJ(1) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaLTE.^2);

if doLSPT
    %% compute the cost for the ALTERED Loss of support scenario
    deltaLoss = 0.01.*LossSupport;
    
    %1) Compute vy_plus
    vy_plus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LossSupport+deltaLoss,LTE);
    vy_plus = vy_plus.*1e3;   %pass from m/sec to mm/sec
    %2) compute SSE_plus
    SSE_plus = vy_plus - TSD_meas;
    SSE_plus = sum(SSE_plus.^2);
    
    %3) compute vy_minus
    vy_minus = getDeflectionVelocity_joint(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LossSupport-deltaLoss,LTE);
    vy_minus = vy_minus.*1e3;   %pass from m/sec to mm/sec
    %4) compute SSE_minus
    SSE_minus = vy_minus - TSD_meas;
    SSE_minus = sum(SSE_minus.^2);
    
    %5) compute gradJ(1) and HessianJ(1)
    gradJ(2) = (SSE_plus - SSE_minus)./(2.*deltaLoss);
    HessianJ(2) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaLoss.^2);
else
    gradJ(2) = 0;
    HessianJ(2) = 0;
end

end %<--- endfunction