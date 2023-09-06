function [E,k,G,LTE,SSE] = backCalc_shortSlab(TSD_points,TSDvy,TSDvx,estimateC,LOAD,pressure,h,nu,subgrade_type,verboseness)
% function = backCalc_shortSlab(TSD_points,TSDvy,TSDvx,estimateC,h,kValue,G,E,nu,LOAD,pressure,verboseness)


%% Front-end function to run a single back-calculation over deflection velocity data for a short slab.
% This back-calc problem solves for c, LTE, and the subgrade loss of support
% with the new implementation based on deflection velocities and using Gradient Descent.
%
%INPUT
%	TSD_points: location of the TSD sensors where the vy values are taken [m]
%   TSDvy:      vertical deflection velocity measurements at the joint's vicinity [m/sec]
%   TSDvx:      TSD traveling speed at the time of measurement [m/sec] - SCALAR!
%   estimateC:  approximate location of the transvese joint from the peak in TSD vy [scalar, m]
%   relative to the TSD wheel [m]
%   h:          concrete slab thickness [m]
%   nu:         concrete's Poisson coef. Dimless, default 0.19-0.23
%   LOAD:       applied LOAD by the TSD 1/2 axle [N]
%   Pressure:   TSD tire pressure [N/m2]
%   subgrade_type: Boolean: 0 = Winkler foundation (G = 0), 1 = Pasternak
%   foundation (G ~=0)
%   verboseness BOOLEAN: 1 = The function reports it's progress all
%   throughout. 0 = silent solver.
%
%OUTPUT
%   kValue:     subgrade's modulus of reaction [N/m3 = 1e9* MPa/mm]
%   G:          subgrade's shear strength [N/m] - for Winkler foundation, G = 0
%   E:          concrete slab elastic modulus [N/m2]
%   LTE:        joint's LTE index. [0-1 value]
%   SSE:        minimum value of SSE achieved for the optimizing [LTE,
%   LossSupport]

%V1.0 - 2023-04-25 Giuietta or not.
% First implementation based on back-calculation solver v2023-01-25.

%% Code begins
TSD_points = TSD_points(:);
TSDvy = TSDvy(:);
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

tryoutLTE = 0.90;

max_iter = 2000;

if verboseness
    fprintf('\t \t Launching short slab back-calculation for k,E,G, and LTE at the joint \n')
end

gradDescentLR = 0.1;
number_of_tryouts = length(tryoutLTE); %<-%this variable is here for QC/QA testing (so that multiple descent attempts can be tried in a single call to this function)

historyE = zeros(max_iter+1,1,number_of_tryouts);
historyk = zeros(max_iter+1,1,number_of_tryouts);
historyG = zeros(max_iter+1,1,number_of_tryouts);
historyLTE = zeros(max_iter+1,1,number_of_tryouts);

%the history gradient's number of columns is the number of decision variables [4] plus 1 (to accomodate L2 size of historyGrad) 
historyGrad = zeros(max_iter,5,number_of_tryouts);
historyHess = zeros(max_iter,4,number_of_tryouts);
historySSE  = zeros(max_iter+1,1,number_of_tryouts);

LossSupport = 1;

for zz = 1:number_of_tryouts
    %% initialize variables
    E = tryoutE(zz);  
    k = tryoutk(zz);
    G = tryoutG(zz);    
    LTE = tryoutLTE(zz);
    
    %% launch the iterative loop that would optimize the values of h,E,g,k (continuous slab ahead of the joint)
    stop = 0;
    iter = 0;
    SSE_old = sum(TSDvy.^2);

    historyLTE(1,1,zz)  = LTE;
    historyE(1,1,zz) = E;
    historyk(1,1,zz) = k;
    historyG(1,1,zz) = G;

    historySSE(1,1,zz)     = SSE_old;

    while (~stop && iter <max_iter)
        iter = iter + 1;       
        %stage 1: compute simulated deflection velocity for given c, LTE, LossSpt.
        W_prima = getDeflectionVelocity_joint(TSD_points,yDomain,TSDvx,h,E,nu,k,G,LOAD,pressure,estimateC,LossSupport,LTE);
        W_prima = W_prima(:);        
        
        %stage 3 - compare iterated deflection slope versus the TSD record
        %use an SSE-based method
        auxSSE = (1e3.*TSDvy - 1e3.*W_prima).^2;
        SSE = sum(auxSSE);        
        
        if iter/100 == floor(iter/100) && verboseness
            fprintf(' \t iteration number %g \n',iter)
            fprintf(' \t joint Loss of support and LTe index are: %g, %g, \n',LossSupport,LTE)
            fprintf(' \t SSE value %g \n',SSE)
        end

        %store the SSE of the current iteration.
        if iter >1
            SSE_old = historySSE(iter);
        end

        %test SSE. If good [SSE stagnated], stop descent; if not, descend for the next
        %iteration.
        if (SSE/SSE_old<1.0001 && SSE/SSE_old>0.9999)
            %satisfactory result, the SSE stagnated! - force the loop to stop
            stop = 1;
        else
            %% do the iterative step
            %Careful!~ here TSD_meas is by default in m/sec (don't need to
            %convert units!)
            [gradJ,Hjj] = computeGradient_shortSlabVy(TSD_points,yDomain,TSDvy,h,nu,E,G,k,estimateC,LTE,LossSupport,LOAD,pressure,TSDvx);
            gradJ = gradJ(:);
            Hjj = Hjj(:);
            historyGrad(iter,1:end-1,zz) = gradJ';
            historyGrad(iter,end,zz) = sqrt(gradJ'*gradJ);
            historyHess(iter,:,zz) = Hjj;
            
            %stability update V2021-12-10: Take the ABS. VALUE of the Hjj
            %entries to prevent the gradient from ascending - case your target function is non-convex [such as SSE when far off, the target, it unfolds to concavity].
            %See Zhang et al., 2021, Pg. 459
            Hjj = abs(Hjj);
            
            %Update k,E,G, LTE based on the descent.             
            k = k - gradDescentLR.*(1./Hjj(1)).*gradJ(1);
            E = E - gradDescentLR.*(1./Hjj(2)).*gradJ(2);            
            if subgrade_type == 0
                G = 0;
            else
                G  = G - gradDescentLR.*(1./Hjj(3))*gradJ(3);
            end
            LTE = LTE - gradDescentLR.*(1./Hjj(4)).*gradJ(4);

            %stability update v2021-12-14
            %if either k or E go to the negative realm, force them back
           if E<0 || E>1e12
                E = 1e10; %if needed, reset E to 1x10^10N/m2 [10 GPa]
            end
            if k<0
                k = 0.033.*1e9; %if needed, reset k to 0.033 MPa/mm [~140 PCI]
            end
            
            %stability control - LTE is BOUNDED BETWEEN 0 AND 1,
            if LTE > 1
                LTE = 1;
            end
            if LTE < 0
                LTE = 0.01;
            end

            %add the record of LTE, Loss, SSE
            historyk(iter+1,1,zz) = k;
            historyE(iter+1,1,zz) = E;     
            historyG(iter+1,1,zz) = G;            
            historyLTE(iter+1,1,zz) = LTE;     
            
            historySSE(iter+1,1,zz) = SSE;  
            
        end %end iterative descent step
    end %end while.
end  %end number_of_tryouts.

%% chop outcomes, remove the zero values of hitsoryLte and historyLoss, SSE Grad
if number_of_tryouts == 1
    historyLTE = historyLTE(historyLTE>0,:,:);
    historySSE = historySSE(historySSE>0,:,:);
end
%% debug zone - plot the domainSSE

% disp('End debug zone')