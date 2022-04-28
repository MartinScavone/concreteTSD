function [gradJ,HessianJ] = computeGradient_JointedVy(xDomain,yDomain,TSD_meas,h,nu,E,G,k,c,LTE,LOAD,pressure,vx)
%function [gradJ,HessianJ] = computeGradient_JointedVy(xDomain,yDomain,TSD_meas,h,nu,E,G,k,c,LTE,LOAD,pressure,vx)
%
%Auxiliary function to the back-calculation tool. 
%Compute the gradient of the cost function for the deflection velocity gradient descent.
%THIS FUNCTION COMPUTES THE GRADIENT FOR THE JOINTED SLAB, WHERE I'M
%BACK-CALCUALTING FOR C, LTE.
%The cost function to minimize is SSE - [TSD_meas - vy(c,LTE|k,E,G)]^2
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
%   c:          distance between load center and joint [m]
%   LTE:        load transfer efficiency index at the joint (deflection-based).
%   LOAD:       amount of load [newtons]
%   pressure:   load pressure [N/m2]
%   vx:         TSD travel speed [scalar, m/sec]  
%
%OUTPUT:
%   gradJ = (approximated) gradient of the cost function for the 6
%   variables, at their current values      
%      gradJ[1]: partialJ/partial_c 
%      gradJ[2]: partialJ/partial_LTE
%     
%       HessianJ[1] = partial2J/partial_c2
%       HessianJ[2] = partial2J/partial_LTE2
%
%release candidate V2022-05-01

%% CODE BEGINS
% pass TSD_meas to mm/sec (because my SSE is all throughout based on vy in
% mm/sec]
TSD_meas = TSD_meas.*1e3;

%Initialize output
gradJ = zeros(2,1);
HessianJ = zeros(2,1);

%% compute base case [needed for HJJ]

%1) compute vy_base
vy_base = getDeflectionVelocity_singleProfile(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LTE);
vy_base = vy_base.*1e3;  %pass from m/sec to mm/sec

%2) compute SSE base
SSE_base = vy_base - TSD_meas;
SSE_base = sum(SSE_base.^2);

%% compute the cost for the ALTERED c scenario
deltaC = 0.01.*c;

%1) Compute vy_plus
vy_plus = getDeflectionVelocity_singleProfile(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c+deltaC,LTE);
vy_plus = vy_plus.*1e3;  %pass from m/sec to mm/sec
%2) compute SSE_plus
SSE_plus = vy_plus - TSD_meas;
SSE_plus = sum(SSE_plus.^2);

%3) compute vy_minus
vy_minus = getDeflectionVelocity_singleProfile(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c-deltaC,LTE);
vy_minus = vy_minus.*1e3;  %pass from m/sec to mm/sec
%4) compute SSE_minus
SSE_minus = vy_minus - TSD_meas;
SSE_minus = sum(SSE_minus.^2);

%5) compute gradJ(1) and HessianJ(1)
gradJ(1) = (SSE_plus - SSE_minus)./(2.*deltaC);
HessianJ(1) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaC.^2);

%% compute the cost for the ALTERED LTE scenario
deltaLTE = 0.01.*LTE;

%1) Compute vy_plus
vy_plus = getDeflectionVelocity_singleProfile(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LTE+deltaLTE);
vy_plus = vy_plus.*1e3;   %pass from m/sec to mm/sec
%2) compute SSE_plus
SSE_plus = vy_plus - TSD_meas;
SSE_plus = sum(SSE_plus.^2);

%3) compute vy_minus
vy_minus = getDeflectionVelocity_singleProfile(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LTE-deltaLTE);
vy_minus = vy_minus.*1e3;   %pass from m/sec to mm/sec
%4) compute SSE_minus
SSE_minus = vy_minus - TSD_meas;
SSE_minus = sum(SSE_minus.^2);

%5) compute gradJ(1) and HessianJ(1)
gradJ(2) = (SSE_plus - SSE_minus)./(2.*deltaLTE);
HessianJ(2) = (SSE_plus - 2.*SSE_base + SSE_minus)./(deltaLTE.^2);

end %<--- endfunction