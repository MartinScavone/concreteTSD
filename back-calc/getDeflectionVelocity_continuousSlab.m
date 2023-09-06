function [vy] = getDeflectionVelocity_continuousSlab(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure)
%function [vy] = getDeflectionVelocity_continuousSlab(xDomain,yDomain,vx,h,E,nu,k,G,LOAD,pressure)
%solve the vertical component of deflection velocity for a single set of
%TSD measurements for the continuous (infinite) slab case.
%In this case, vy = vx.*Slope!
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
%   vx:         TSD travel speed [scalar, m/sec]
%OUTPUT
%   vy = matrix with the vertical velocity component [m/sec]
%at the moment of survey for each location xDomain,yDomain
%
%V 1.0 - 2022-04-21 Maconha Alleman Aftermath
%   Built from old deflectionVelocity v0413 function.
%   Compute vy using finite differences centered at the measurement
%   location. Use a DELTA displacement of 0.01m
%
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
wTplus = computeContinuousSlabDeflectionBasin(xDomain-deltaX,yDomain,h,E,nu,k,G,LOAD,pressure);
wTminus = computeContinuousSlabDeflectionBasin(xDomain+deltaX,yDomain,h,E,nu,k,G,LOAD,pressure);

%solve vy as centered difference - hence the 2*deltaT denominator!
vy = (wTplus-wTminus)./(2.*deltaT);      

end