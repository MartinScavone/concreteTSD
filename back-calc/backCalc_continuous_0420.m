function [E,k,G,historyE,historyk,historyG,historySSE] = backCalc_continuous_0420(TSD_points,TSD_meas,TSDvx,Load,pressure,h,nu,subgrade_type,verboseness)
% function [E,k,G,historyE,historyk,historyG,historySSE] = backCalc_continuous_0420(TSD_points,TSD_meas,TSDvx,Load,pressure,h,nu,subgrade_type,verboseness)
%% Front-end function to run a single back-calculation over deflection velocity data for the continuous slab w/o joint.
% This back-calc problem solves for E, k G; with the new implementation
% based on deflection velocities.
%
%INPUT
%	TSD_points: locations (w.r.t the TSD rear wheel) where the slope defl. measurements are taken [m]
%   TSD_meas: 	deflection velocity measurements from the TSD [m/sec]
%   TSDvx:      TSD travel speed [m/sec]
%   Load: 		half-axle load [Newtons]
%   Pressure: 	half-axle tire pressure [Pa]
%   h: 			concrete slab thickness [m]
%   nu: 		concrete Poisson coefficient [default 0.20]
%   subgrade_type: Boolean: 0 = Winkler foundation (G = 0), 1 = Pasternak
%   foundation (G ~=0)
%   verboseness BOOLEAN: 1 = The function reports it's progress all
%   throughout. 0 = silent solver.
%
%OUTPUT
%   E: Back-calculated slab's Young modulus [N/m2]
%   k: Back-calculated subgrade's Mod. of reaction [N/m3]
%   G: subgrade's shear modulus [N/m2] (If subgrade_type = 0), return is G = 0
%   historyE: evolution of back-calc E values
%   historyK: evolution of back-calc k values
%   historyG: evolution of back-calc G values
%   historySSE: Evolution of target function [SSE between measurements and back-calc defl. basins]
%
%release candidate v2022-05-01


%% CALCULATE DEFLECTIONS
TSD_points = TSD_points(:);
TSD_meas = TSD_meas(:);
yDomain = 0;

%% LAUNCH BACK-CALCULATION!
tryoutE = 25e9; %Young modulus for concrete [N/m2]
tryoutk = 100; %defaults in PCI
tryoutk = tryoutk./3.684.*1e6;  %same, now in N/m3 [thanks, Google]

if subgrade_type == 0
        tryoutG = 0;  %G = 0 for Winkler foundation, otherwise, it's Pasternak foundation.
else
        tryoutG = 6e8;
end
number_of_tryouts = length(tryoutE);
max_iter = 2000;

if verboseness
    fprintf('\t \t Launching back-calculation for k, E, G \n')
end

gradDescentLR = 0.1;

historyE = zeros(max_iter+1,1,number_of_tryouts);
historyk = zeros(max_iter+1,1,number_of_tryouts);
historyG = zeros(max_iter+1,1,number_of_tryouts);

historyGrad = zeros(max_iter,4,number_of_tryouts);
historyHess = zeros(max_iter,3,number_of_tryouts);
historySSE = zeros(max_iter+1,1,number_of_tryouts);

for zz = 1:number_of_tryouts
    %% initialize variables
    E = tryoutE(zz);  
    k = tryoutk(zz);
    G = tryoutG(zz);    
    
    %% launch the iterative loop that would optimize the values of h,E,g,k (continuous slab ahead of the joint)
    stop = 0;
    iter = 0;
    SSE_old = sum(TSD_meas.^2);

    historyE(1,1,zz) = E;
    historyk(1,1,zz) = k;
    historyG(1,1,zz) = G;
    historySSE(1) = SSE_old;

    while (~stop && iter <max_iter)
        iter = iter + 1;       
        %stage 1: compute simulated deflection veloctiy\
        W_prima = getDeflectionVelocity_continuousSlab(TSD_points,yDomain,TSDvx,h,E,nu,k,G,Load,pressure);
        W_prima = W_prima(:);        
        
        %stage 3 - compare iterated deflection slope versus the TSD record
        %use an SSE-based method
        auxSSE = (1e3.*TSD_meas - 1e3.*W_prima).^2;
        SSE = sum(auxSSE);        
        
        if iter/100 == floor(iter/100) && verboseness
            fprintf(' \t iteration number %g \n',iter)
            fprintf(' \t E, k, G values are: %g, %g, %g \n',E,k,G)
            fprintf(' \t SSE value %g \n',SSE)
        end

        %stage 4 - use the results from stage 4 to improve the prediction
        %%%GRADIENT DESCENT -> V0.12: do Grad Descent over l, D, g!!!!
        if iter >1
            SSE_old = historySSE(iter);
        end

        if (SSE/SSE_old<1.0001 && SSE/SSE_old>0.9999)
            %satisfactory result, the SSE stagnated! - force the loop to stop
            stop = 1;
        else
            %% do the iterative step
            %Careful!~ here TSD_meas is by default in m/sec (don't need to
            %convert units!)
            [gradJ,Hjj] = computeGradient_ContinuousVy(TSD_points,yDomain,TSD_meas,h,nu,E,G,k,Load,pressure,TSDvx);
            gradJ = gradJ(:);
            Hjj = Hjj(:);
            historyGrad(iter,1:3,zz) = gradJ';
            historyGrad(iter,4,zz) = sqrt(gradJ'*gradJ);
            historyHess(iter,:,zz) = Hjj;
            
            %stability update V2021-12-10: Take the ABS. VALUE of the Hjj
            %entries to prevent the gradient from ascending - case your target function is non-convex [such as SSE when far off, the target, it unfolds to concavity].
            %See Zhang et al., 2021, Pg. 459
            Hjj = abs(Hjj);
            
            %Update k, E, G based on the descent.            
            k = k - gradDescentLR.*(1./Hjj(1)).*gradJ(1);
            E = E - gradDescentLR.*(1./Hjj(2)).*gradJ(2);
            
            %stabilty update V2021-08-23. If G = 0 since the beginning, don't update g cause it causes a crash! 
            %[I can control this by checking on subgrade_type]
            if subgrade_type == 0
                G = 0;
            else
                G  = G - gradDescentLR.*(1./Hjj(3))*gradJ(3);
            end
            
            %stability update v2021-12-14
            %if either k or E go to the negative realm, force them back
           if E<0 || E>1e12
                E = 1e10; %if needed, reset E to 1x10^10N/m2 [10 GPa]
            end
            if k<0
                k = 0.033.*1e9; %if needed, reset k to 0.033 MPa/mm [~140 PCI]
            end
            historyk(iter+1,1,zz) = k;
            historyE(iter+1,1,zz) = E;     
            historyG(iter+1,1,zz) = G;
            historySSE(iter+1,1,zz) = SSE;  
            
            %stop = 0
        end %end iterative descent step
    end %end while.
end  %end number_of_tryouts.

%% chop outcomes, remove the zero values of historyk, E, G, SSE Grad
historyk = historyk(historyk>0,:);
historyE = historyE(historyE>0,:);
if subgrade_type ~= 0
    historyG = historyG(historyG>0,:);
else
    historyG = historyG(historyE>0,:);
end    
historySSE = historySSE(historySSE>0,:);
% historyGrad = historyGrad(historyGrad(:,end)>0,:);
% historyHess = historyHess(historyH(:,end)>0,:);

