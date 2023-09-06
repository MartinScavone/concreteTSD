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
%   G: subgrade's shear modulus [N/m] (If subgrade_type = 0), return is G = 0
%   historyE: evolution of back-calc E values
%   historyK: evolution of back-calc k values
%   historyG: evolution of back-calc G values
%   historySSE: Evolution of target function [SSE between measurements and back-calc defl. basins]
%
%V 3.0 - 2022-04-21 Maconha Alleman Aftermath
%   Code overhauled for the new implmentation based entirely on vy data.
%   Code simplified and dependencies streamlined.

%V 2.2 - 2022-04-11
%   Updated the tire width from 0.68 to 0.50m (width of a heavy truck
%   dual-wheel single axle -as measured)
%V 2.1.4 - March Negative 2022-03-02
%   Remove the correction from v2.1.3 (it will crash a test on simulated
%   [with correct sign] TSD data!
%   Prompt the end-user to instead add the -1 multiplier when passing real
%   TSD data into the back-calc engine.
%V 2.1.3 - Negative - 2021-12-10
%   Multiply the TSD slope-defl measurements by -1 to match convention with
%   the back-calc responses.
%   Re-written to match the testScript v0.3 (back-calc directly on k, E;
%   stabilized by Hessian term in descent)%
%V 2.0 - Virgin's day - 2021-12-09
% 	Compatibility with v2.0 back-calc front-end
%	Added the TSD_point locations as input to the back-calc engine (to allow for more flexible input).	%
%V 1.1 - MnRoad - 2021-12-06
%   Removed a "clear all" statement at the beginning, which killed the
%   input variables.%
%V 1.0 Nostalgia Night (2021-08-24)
%   First functional version. Based on code v0.2

%% preparation


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
number_of_tryouts = length(tryoutE); %<-%this variable is here for QC/QA testing (so that multiple descent attempts can be tried in a single call to this function)

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

%% debug zone - plot the domainSSE
%         % graphical output
%         startE = historyE(1);
%         startK = historyk(1);
%         
%         figure(4)
%         subplot(2,2,1)
%         plot(historyE(historyE>0),historyk(historyE>0))
%         hold on
%         plot(startE,startK,'r*','markersize',10)
%         hold off
%         xlabel('E [N/m2]')
%         ylabel('k [N/m3]')
%         grid on
% 
%         subplot(2,2,2)
%         nSE = sum(historySSE>0);
%         plot(1:1:nSE,historySSE(historySSE>0),'-k+','markersize',8)
%         xlabel('epoch')
%         ylabel('SSE')
% 
%         subplot(2,2,3)
%         nSE = sum(historyGrad(:,4)>0);
%         plot(1:1:nSE,historyGrad(1:nSE,4),'-k+','markersize',8)
%         xlabel('epoch')
%         ylabel('L2 norm of gradient of SSE')
% 
%         subplot(2,2,4)
%         nSE = sum(historySSE>0);
%         yyaxis left
%         plot(1:nSE,historyE(1:nSE),'b')
%         grid
%         xlabel('iteration')
%         ylabel('E [N/m2]')
% 
%         yyaxis right
%         plot(1:nSE,historyk(1:nSE),'r')
%         ylabel('k [N/m3]')
% 
%         %% comparison plot - SSE(E,k) 
%         % domainE = 0.25e10:0.25e10:3.5e10;
%         domainE = 1e10:0.25e10:6e10;
%         domainK = 1e7:0.25e7:10e7;
%         domainSSE = zeros(length(domainE),length(domainK));
%         disp('computing domain for SSE')
% 
%         for i = 1:length(domainE)
%             for j = 1:length(domainK)
%                 l = ((domainE(i)*h^3)/(12*(1-nu^2)*domainK(j)))^0.25;
%                 DD= domainK(j).*l^4;
%                 g = G*l^2/(2*DD);
%                 W_prima = compute_w_dot(TSD_points,yDomain,sDomain,pressure,a,b,domainK(j),g,l);
%                 W_prima = W_prima(:,yDomain==0);  
%                 %updated v2021-12-10
%                 W_prima = W_prima.*1e6;        
%                 domainSSE(i,j) = sum((TSD_meas - W_prima).^2);
%             end
%         end
% 
% 
%         %%
%         figure(2)
%         subplot(1,2,1)
%         surf(domainK, domainE,real(domainSSE))
%         xlabel('k [N/m3]')
%         ylabel('E [N/m2]')
%         zlabel ('SSE')
%         hold on
%         plot(startK,startE,'r*','markersize',15)
%         plot3(historyk(historyE>0),historyE(historyE>0),historySSE(historyE>0),'y','linewidth',2)
%         xlim([0,8e7]);
%         ylim([0,3.5e10]);
%         hold off
% 
%         subplot(1,2,2)
%         contour(domainK, domainE,real(domainSSE))
%         xlabel('k [N/m3]')
%         ylabel('E [N/m2]')
%         zlabel ('SSE')
%         grid on
%         hold on
%         plot(startK,startE,'r*','markersize',15)
%         plot(historyk(historyE>0),historyE(historyE>0),'k','linewidth',2)
%         plot(historyk(end),historyE(end),'k*','markersize',15)        
%         hold off
%         colorbar('southoutside')
%         legend('SSE','start','descent','target')