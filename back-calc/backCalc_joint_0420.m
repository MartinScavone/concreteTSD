function [LTE,LossSupport,SSE,historyLTE,historyLossSpt,historySSE] = backCalc_joint_0420(TSD_points,TSDvy,TSDvx,estimateC,h,kValue,G,E,nu,LOAD,pressure,doLSPT,verboseness)
% function [LTE,LossSupport,SSE,historyLTE,historyLossSpt,historySSE] = backCalc_joint_0420(TSD_points,TSDvy,TSDvx,estimateC,h,kValue,G,E,nu,LOAD,pressure,doLSPT,verboseness)
%% Front-end function to run a single back-calculation over deflection velocity data for the jointed slab.
% This back-calc problem solves for LTE, and the subgrade loss of support
% with the new implementation based on deflection velocities and using Gradient Descent.
%
%INPUT
%	TSD_points: location of the TSD sensors where the vy values are taken [m]
%   TSDvy:      vertical deflection velocity measurements at the joint's vicinity [m/sec]
%   TSDvx:      TSD traveling speed at the time of measurement [m/sec] - SCALAR!
%   estimateC:  approximate location of the transvese joint from the peak  in TSD vy [scalar, m] relative to the TSD wheel [m]
%   h:          concrete slab thickness [m]
%   kValue:     subgrade's modulus of reaction [N/m3 = 1e9* MPa/mm]
%   G:          subgrade's shear strength [N/m] - for Winkler foundation, G = 0
%   E:          concrete slab elastic modulus [N/m2]
%   nu:         concrete's Poisson coef. Dimless, default 0.19-0.23
%   LOAD:       applied LOAD by the TSD 1/2 axle [N]
%   Pressure:   TSD tire pressure [N/m2]
%   doLSPT:     Boolean: 1 to solver for loss of support at the subgrade, 0 if not.
%   verboseness: Boolean: 1 if you want the solver to be verbose, 0 otherwise
%
%OUTPUT
%   LTE:        joint's LTE index. [0-1 value]
%   LossSupport:Subgrade's loss of support. 0-1 value]
%   SSE:        minimum value of SSE achieved for the optimizing [LTE,
%   LossSupport]
%   hisstoryLTE:all the tried values for the joint's LTE index.    
%   historyLoss:all the tried values for the Loss of support
%   historySSE: Error domain SSE(domainLTE,domainLoss)

%V2.1 - 2023-0405 Spring
% Add a boolean input variable "doLSPT" to control whether or not to
% compute Loss of support (and roll back to original back-calc formulation,
% w/o LSPT)
%V2.0 - 2023-01-25 Flu.
% Experimental version to solve the jointed slab problem with loss of support of the
% subgrade [Loss] and joint's LTE index by Gradient Descent.


%% Code begins
TSD_points = TSD_points(:);
TSDvy = TSDvy(:);
yDomain = 0;

%% LAUNCH BACK-CALCULATION!
% tryoutLTE = [0.50 0.90];
% tryoutLSPT= [0.40 0.80];

if doLSPT
    tryoutLTE = 0.90;
    tryoutLSPT= 0.80;
else
    tryoutLTE = 0.90;
    tryoutLSPT = 1;  %no LSPT calculation, assume subgrade with full support.
end
max_iter = 2000;

if verboseness
    fprintf('\t \t Launching back-calculation for LTE and loss of support at the joint \n')
end

gradDescentLR = 0.1;
number_of_tryouts = length(tryoutLTE); %<-%this variable is here for QC/QA testing (so that multiple descent attempts can be tried in a single call to this function)

historyLTE = zeros(max_iter+1,1,number_of_tryouts);
historyLossSpt = zeros(max_iter+1,1,number_of_tryouts);

historyGrad = zeros(max_iter,3,number_of_tryouts);
historyHess = zeros(max_iter,2,number_of_tryouts);
historySSE  = zeros(max_iter+1,1,number_of_tryouts);

for zz = 1:number_of_tryouts
    %% initialize variables
    LTE = tryoutLTE(zz);
    LossSupport = tryoutLSPT(zz);
    
    %% launch the iterative loop that would optimize the values of LTE and LSPT
    stop = 0;
    iter = 0;
    SSE_old = sum(TSDvy.^2);

    historyLTE(1,1,zz)  = LTE;
    historyLossSpt(1,1,zz) = LossSupport;
    historySSE(1,1,zz)     = SSE_old;

    while (~stop && iter <max_iter)
        iter = iter + 1;       
        %stage 1: compute simulated deflection velocity for given c, LTE, LossSpt.
        W_prima = getDeflectionVelocity_joint(TSD_points,yDomain,TSDvx,h,E,nu,kValue,G,LOAD,pressure,estimateC,LossSupport,LTE);
        W_prima = W_prima(:);        
        
        %stage 3 - compare iterated deflection slope versus the TSD record
        %use an SSE-based method
        auxSSE = (1e3.*TSDvy - 1e3.*W_prima).^2;
        SSE = sum(auxSSE);        
        
        if iter/100 == floor(iter/100) && verboseness
            fprintf(' \t iteration number %g \n',iter)
            fprintf(' \t joint Loss of support and LTE index are: %g, %g, \n',LossSupport,LTE)
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
            %--update v2023-04-05 -> Pass doLSPT here too!
            [gradJ,Hjj] = computeGradient_JointedVy(TSD_points,yDomain,TSDvy,h,nu,E,G,kValue,estimateC,LTE,LossSupport,LOAD,pressure,TSDvx,doLSPT);
            gradJ = gradJ(:);
            Hjj = Hjj(:);
            historyGrad(iter,1:end-1,zz) = gradJ';
            historyGrad(iter,end,zz) = sqrt(gradJ'*gradJ);
            historyHess(iter,:,zz) = Hjj;
            
            %stability update V2021-12-10: Take the ABS. VALUE of the Hjj
            %entries to prevent the gradient from ascending - case your target function is non-convex [such as SSE when far off the target, it unfolds to concavity].
            %See Zhang et al., 2021, Pg. 459
            Hjj = abs(Hjj);
            
            %Update c, LTE, LossSupport based on the descent. 
            if doLSPT
                LTE = LTE - gradDescentLR.*(1./Hjj(1)).*gradJ(1);
                LossSupport = LossSupport - gradDescentLR.*(1./Hjj(2)).*gradJ(2);
            else
                LTE = LTE - gradDescentLR.*(1./Hjj(1)).*gradJ(1);
%                 LossSupport = LossSupport;
            end

            %stability control - LTE AND LOSS ARE BOUNDED BETWEEN 0 AND 1,
            if LTE > 1
                LTE = 1;
            end
            if LTE < 0
                LTE = 0.01;
            end

            if LossSupport<0
                LossSupport = 0.02;
            end
            if LossSupport > 1
                LossSupport = 1;
            end

            %add the record of LTE, Loss, SSE
            historyLTE(iter+1,1,zz) = LTE;     
            historyLossSpt(iter+1,1,zz) = LossSupport;
            historySSE(iter+1,1,zz) = SSE;  
            
        end %end iterative descent step
    end %end while.
end  %end number_of_tryouts.

%% chop outcomes, remove the zero values of hitsoryLte and historyLoss, SSE Grad
if number_of_tryouts == 1
    historyLTE = historyLTE(historyLTE>0,:,:);
    historyLossSpt = historyLossSpt(historyLossSpt>0,:,:);
    historySSE = historySSE(historySSE>0,:,:);
end
%% debug zone - plot the domainSSE
% % graphical output
% startLTE = historyLTE(1,:,1);
% startLoss = historyLossSpt(1,:,1);
% startLTE2 = historyLTE(1,:,2);
% startLoss2 = historyLossSpt(1,:,2);
% 
% figure(4)
% subplot(2,2,1)
% plot(historyLTE(historyLTE(:,:,1)>0,:,1),historyLossSpt(historyLTE(:,:,1)>0,:,1))
% hold on
% plot(startLTE,startLoss,'r*','markersize',10)
% hold off
% xlabel('LTE')
% ylabel('Loss of support')
% grid on
% 
% subplot(2,2,2)
% nSE = sum(historySSE(:,:,1)>0);
% plot(1:1:nSE,historySSE(historySSE(:,:,1)>0),'-k+','markersize',8)
% xlabel('epoch')
% ylabel('SSE')
% 
% subplot(2,2,3)
% nSE = sum(historyGrad(:,end,1)>0);
% plot(1:1:nSE,historyGrad(1:nSE,end),'-k+','markersize',8)
% xlabel('epoch')
% ylabel('L2 norm of gradient of SSE')
% 
% subplot(2,2,4)
% yyaxis left
% plot(1:nSE,historyLTE(1:nSE),'b')
% grid
% xlabel('iteration')
% ylabel('LTE')
% 
% yyaxis right
% plot(1:nSE,historyLossSpt(1:nSE),'r')
% ylabel('subgrade loss of support')
%
%% plot the domain of LTE and Loss of support
% domainLTE = 0.10:0.02:1;
% domainLoss = 0.20:0.05:1;
% domainSSE = zeros(length(domainLTE),length(domainLoss));
% disp('computing domain for SSE')
% 
% for i = 1:length(domainLTE)
%   for j = 1:length(domainLoss)          
%       W_prima = getDeflectionVelocity_joint(TSD_points,yDomain,TSDvx,h,E,nu,kValue,G,LOAD,pressure,estimateC,domainLoss(j),domainLTE(i));
%       W_prima = W_prima(:); 
%       
%       auxSSE = (1e3.*TSDvy - 1e3.*W_prima).^2;
%       domainSSE(i,j) = sum(auxSSE);        
%   end
% end
% 
% figure(2)
% subplot(1,2,1)
% surf(domainLoss, domainLTE,real(domainSSE))
% set(gca,'fontsize',16)
% xlabel('Loss of support')
% ylabel('LTE')
% zlabel ('SSE')
% hold on
% plot(startLoss,startLTE,'r*','markersize',15)
% % historyLTE(historyLTE(:,:,1)>0,:,1),historyLossSpt(historyLTE(:,:,1)>0
% plot3(historyLossSpt(historyLTE(:,:,1)>0,:,1),historyLTE(historyLTE(:,:,1)>0,:,1),historySSE(historyLTE(:,:,1)>0,:,1),'y','linewidth',2)
% plot(startLoss2,startLTE2,'r*','markersize',15)
% plot3(historyLossSpt(historyLTE(:,:,2)>0,:,2),historyLTE(historyLTE(:,:,2)>0,:,2),historySSE(historyLTE(:,:,2)>0,:,2),'g','linewidth',2)
% hold off
% 
% subplot(1,2,2)
% contour(domainLoss, domainLTE,real(domainSSE))
% set(gca,'fontsize',16)
% xlabel('Loss of support')
% ylabel('LTE')
% zlabel ('SSE')
% grid on
% hold on
% plot(startLoss,startLTE,'r*','markersize',15)
% plot(historyLossSpt(historyLTE(:,:,1)>0,:,1),historyLTE(historyLTE(:,:,1)>0,:,1),'k','linewidth',2)
% plot(startLoss2,startLTE2,'r*','markersize',15)
% plot(historyLossSpt(historyLTE(:,:,2)>0,:,2),historyLTE(historyLTE(:,:,2)>0,:,2),'b','linewidth',2)
% groundTruthLTE = 0.75;
% groundTruthLSPT = 0.60;   
% plot(groundTruthLSPT,groundTruthLTE,'k*','markersize',15)        
% hold off
% colorbar('southoutside')
% legend('SSE','start_1','descent_1','start_2','descent_2','target')  
% 
% disp('End debug zone')