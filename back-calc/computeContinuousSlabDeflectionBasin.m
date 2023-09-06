function w = computeContinuousSlabDeflectionBasin(xDomain,yDomain,h,E,nu,k,G,LOAD,pressure)
% function w = computeContinuousSlabDeflectionBasin(xDomain,yDomain,h,E,nu,k,G,LOAD,pressure)
%Front-end function to solve the deflection basin of a continuous slab as per
%Van Cauwelaert (2004). This function automates what I had so far been doing on scripts.
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
%OUTPUT
%   w:          deflection basin for each xDomain,yDomain position [m]
%
%V 2.0 - 2023-01-09 TRB
%   Added the Loss of Support estimation (multiplier to k-value) to
%   calculations.
%
%V1.0 2022-04-20 - Maconha Alleman
%First version of this function. Made from scripts v 2022-04-20 and 

%% preparation - calculate loaded area. Assume default 2b = 0.47m
b = 0.47;  %typical width of a 1/2 axle (dual tire) [m]
a = LOAD./(b.*pressure);  %length of distributed load zone [m]
b = b/2;
a = a/2;

%% compute Van C. parameters
LL = ((E.*h.^3)./(12.*(1-nu^2).*k))^0.25;
DD = (E.*h.^3)./(12.*(1-nu^2));
gg = (G.*LL.^2)./(2*DD); 

%% dummy variable for w(x,y) integration
sDomain = 0:0.2:3;
sDomain(1) = 1e-4;

%% This is the normal case (infinite slab, no boundary conditions)
w = compute_w(xDomain,yDomain,sDomain,pressure,a,b,k,gg,LL);

