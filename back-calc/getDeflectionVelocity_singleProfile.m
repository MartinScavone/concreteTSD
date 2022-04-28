function [vy] = getDeflectionVelocity_singleProfile(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LTE)
%function [vy] = getDeflectionVelocity_singleProfile(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure,c,LTE)
%solve the vertical component of deflection velocity for a single set of TSD measurements
%
%INPUT
%   xDomain:    caculation domain, longitudinal direction [m]
%   yDomain:    calculation domain, transverse direction [m]
%   h:          slab thickness [m]
%   E:          slab elastic modulus [N/m2]
%   nu:         slab Poisson coefficient [dimless]
%   k:          subgrade's modulus of reaction [N/m3]
%   G:          subgrade's shear modulus [N/m2]. For Winkler foundation, G =0
%   LOAD:       amount of load [newtons]
%   pressure:   load pressure [N/m2]
%   c:          distance between the joint and the load center [m]
%   LTE:        joint's load transfer efficiency index [dimless]
%   vx:         TSD travel speed [scalar, m/sec]
%OUTPUT
%   vy = matrix with the vertical velocity component [m/sec]
%at the moment of survey for each location xDomain,yDomain
%
%release candidateV2022-05-01
%%

deltaX = 0.01;
deltaT = deltaX./vx;

%initialize vy
vy = zeros(length(xDomain),length(yDomain));

%Solve vy as a finite difference centered at the measurement point.
%Assume that the load is approaching the joint as time passes. Thus the
%xDomain relative to the load reduces by an amount deltaX for the "time
%plus" measurement and increases by deltaX for the "time minus"
%measurement.
%Also, the distance c changes by -/+ deltaX as well (because the load got
%closer to the joint!)

% call auxiliar function to solve the deflection basins
wTplus = computeJointedSlabDeflectionBasin(xDomain-deltaX,yDomain,h,E,nu,k,G,LOAD,pressure,c-deltaX,LTE);
wTminus = computeJointedSlabDeflectionBasin(xDomain+deltaX,yDomain,h,E,nu,k,G,LOAD,pressure,c+deltaX,LTE);

vy = (wTplus-wTminus)./(2.*deltaT);      

end