function w = computeJointedSlabDeflectionBasin(xDomain,yDomain,h,E,nu,k,G,LOAD,pressure,c,LossSupport,LTE)
% function w = computeJointedSlabDeflectionBasin(xDomain,yDomain,h,E,nu,k,G,LOAD,pressure,c,LossSupport,LTE)
%Front-end function to solve the deflection basin of a jointed slab as per
%Van Cauwelaert (2004). This function automates what I had so far been
%doing on scripts.
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
%   LossSupport Loss of support of the subgrade [dimless]
%   LTE:        joint's load transfer efficiency index [dimless]
%
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

%% update v2023-01-09 - add the subgrde loss of support term to k-value (multiply k by the loss-of-support to drop the value)
k = k.*LossSupport;

%% compute Van C. parameters
LL = ((E.*h.^3)./(12.*(1-nu^2).*k))^0.25;
DD = (E.*h.^3)./(12.*(1-nu^2));
gg = (G.*LL.^2)./(2*DD); 

%% dummy variable for w(x,y) integration
sDomain = 0:0.2:3;
sDomain(1) = 1e-4;

ABCDverboseness = 0;%   
w = zeros(length(xDomain),length(yDomain));

%% solve deflection basin
%Apply update v 2022-03-19: Superposition principle
%Check here if the load is stepping on both sides of the joint or not.
%If so, divide the problem into two.

checkPlus  = c>0 && abs(c)<a;   %load on the approach slab but touching the joint
checkMinus = c<0 && abs(c)<a;   %load on the leave slab but touching the joint
    
if checkPlus
    %Load on the approach slab and invading the joint
    %Divide the problem into two chunks of load touching the joint.
    %% problem 1: portion of the load on the approach slab
    a1 = (a + abs(c))./2;
    c1 = a1;
    [AS,BS,CS,DS] = solveABCD_experimental(yDomain,sDomain,pressure,a1,b,k,gg,LL,nu,LTE,c1,ABCDverboseness);
    W1 = zeros(length(xDomain),length(yDomain));
    x1 = find(xDomain<=c1);
    x2 = find(xDomain>c1);
    [Wa,Wb] = compute_wawb(xDomain(x1),yDomain,sDomain,pressure,b,k,gg,LL,AS,BS);
    [Wc,Wd] = compute_wcwd(xDomain(x2),yDomain,sDomain,pressure,b,k,gg,LL,CS,DS);    
    Wx1 = compute_w(xDomain(x1),yDomain,sDomain,pressure,a1,b,k,gg,LL);
    
    W1(x1,:) = Wx1 + Wa + Wb;   %<<--- These deflections are in METERS!!!
    W1(x2,:) = Wc + Wd;         %<<--- Van C. (2004) formulation.      

    %% problem 2: portion of the load on the leave slab
    a2 = (a - abs(c))./2;
    c2 = -a2;%made negative on purpose to highlight that the second chunk of load is a 'reversed' problem [it's on the leave slab]
    %but c2 must pass as a positive value to the ABCD solver!
    [AS,BS,CS,DS] = solveABCD_experimental(yDomain,sDomain,pressure,a2,b,k,gg,LL,nu,LTE,-c2,ABCDverboseness);
    W2 = zeros(length(xDomain),length(yDomain));
    t = -xDomain;
    x1 = find(t<=-c2);
    x2 = find(t>-c2);
    [Wa,Wb] = compute_wawb(t(x1),yDomain,sDomain,pressure,b,k,gg,LL,AS,BS);
    [Wc,Wd] = compute_wcwd(t(x2),yDomain,sDomain,pressure,b,k,gg,LL,CS,DS);    
    Wx1 = compute_w(xDomain(x1),yDomain,sDomain,pressure,a2,b,k,gg,LL);
    
    W2(x1,:) = Wx1 + Wa + Wb;  %<<---  these deflections are in METERS!!!
    W2(x2,:) = Wc + Wd;  %<<--- Van C. (2004) formulation. 
    
    %sum both problems to get the final defl. basin [function output]
    w = W1 + W2;
    
elseif checkMinus
    %Load on the leave slab and invading the joint
    %Divide the problem into two chunks of load touching the joint.
    %% problem 1: portion of the load on the approach slab
    a1 = (a - abs(c))./2;
    c1 = a1;
    [AS,BS,CS,DS] = solveABCD_experimental(yDomain,sDomain,pressure,a1,b,k,gg,LL,nu,LTE,c1,ABCDverboseness);
    W1 = zeros(length(xDomain),length(yDomain));
    x1 = find(xDomain<=c1);
    x2 = find(xDomain>c1);
    [Wa,Wb] = compute_wawb(xDomain(x1),yDomain,sDomain,pressure,b,k,gg,LL,AS,BS);
    [Wc,Wd] = compute_wcwd(xDomain(x2),yDomain,sDomain,pressure,b,k,gg,LL,CS,DS);    
    Wx1 = compute_w(xDomain(x1),yDomain,sDomain,pressure,a1,b,k,gg,LL);
    
    W1(x1,:) = Wx1 + Wa + Wb;   %<<--- These deflections are in METERS!!!
    W1(x2,:) = Wc + Wd;         %<<--- Van C. (2004) formulation.      

    %% problem 2: portion of the load on the leave slab
    a2 = (a + abs(c))./2;
    c2 = -a2;%made negative on purpose to highlight that the second chunk of load is a 'reversed' problem [it's on the leave slab]
    %but c2 must pass as a positive value to the ABCD solver!
    [AS,BS,CS,DS] = solveABCD_experimental(yDomain,sDomain,pressure,a2,b,k,gg,LL,nu,LTE,-c2,ABCDverboseness);
    W2 = zeros(length(xDomain),length(yDomain));
    t = -xDomain;
    x1 = find(t<=-c2);
    x2 = find(t>-c2);
    [Wa,Wb] = compute_wawb(t(x1),yDomain,sDomain,pressure,b,k,gg,LL,AS,BS);
    [Wc,Wd] = compute_wcwd(t(x2),yDomain,sDomain,pressure,b,k,gg,LL,CS,DS);    
    Wx1 = compute_w(xDomain(x1),yDomain,sDomain,pressure,a2,b,k,gg,LL);
    
    W2(x1,:) = Wx1 + Wa + Wb;  %<<---  these deflections are in METERS!!!
    W2(x2,:) = Wc + Wd;  %<<--- Van C. (2004) formulation. 
    
    %sum both problems to get the final defl. basin [function output]
    w = W1 + W2;    
    
else
    %% This is the normal case (load away from the joint), approach slab.
    if c >=0
        [AS,BS,CS,DS] = solveABCD_experimental(yDomain,sDomain,pressure,a,b,k,gg,LL,nu,LTE,c,ABCDverboseness);
        %stage 2: compute deflection
%         W = zeros(length(xDomain),length(yDomain));
        x1 = find(xDomain<=c);
        x2 = find(xDomain>c);
        [Wa,Wb] = compute_wawb(xDomain(x1),yDomain,sDomain,pressure,b,k,gg,LL,AS,BS);
        [Wc,Wd] = compute_wcwd(xDomain(x2),yDomain,sDomain,pressure,b,k,gg,LL,CS,DS);
        Wx1 = compute_w(xDomain(x1),yDomain,sDomain,pressure,a,b,k,gg,LL);
        
        w(x1,:) = Wx1 + Wa + Wb; %<<---  these deflections are in METERS!!!
        w(x2,:) = Wc + Wd;       %<<--- Van C. (2004) formulation. 
%         w = W;
    else
        %% This is the normal case (load away from the joint), leave slab.
        %REVERTED PROBLEM: Revert the boundary conditions
        [AS,BS,CS,DS] = solveABCD_experimental(yDomain,sDomain,pressure,a,b,k,gg,LL,nu,LTE,-c,ABCDverboseness);
        %stage 2: compute deflection 
%         W = zeros(length(xDomain),length(yDomain));
        % AFTER THE JOINT IT'S THE REVERTED PROBLEM, define auxiliary variable t = -x
        t = -xDomain;
        x1 = find(t<=-c);
        x2 = find(t>-c);
        %
        [Wa,Wb] = compute_wawb(t(x1),yDomain,sDomain,pressure,b,k,gg,LL,AS,BS);
        [Wc,Wd] = compute_wcwd(t(x2),yDomain,sDomain,pressure,b,k,gg,LL,CS,DS);        
        Wx1 = compute_w(xDomain(x1),yDomain,sDomain,pressure,a,b,k,gg,LL);
        
        w(x1,:) = Wx1 + Wa + Wb;%<<---  these deflections are in METERS!!!
        w(x2,:) = Wc + Wd;      %<<--- Van C. (2004) formulation.         
        
    end
end
    