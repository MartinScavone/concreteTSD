%% Concrete back-calculation - FRONT-END SCRIPT FOR TSD BACK-CALCULATION
%%- TSD RUN 5cm from MnROAD section #239

%V6.0 - 2023-04-25 - 
% -- Updated front=end to solve the case of pavements with short slabs
% (back-calculation that simultaneously solves the kEG and the LTE of the
% joint)
% -- added numerical solver for LTE (removed correction for 'c'
% -- Experimental: added estimate of subgrae loss of support (LSPT). User must turn variable doLSPT = 1 to do the calculation.
%
%Candidate release updated 2023-09-06

%% STAGE 0 DATA LOAD
%LOAD HERE THE TSD RECORDS - or the BP denoised slope-deflection records
%TSD De-noising procedure front-end code
%
%Copy the bpd + rwl1 front-end code
tic
clc
restoredefaultpath
clear variables
close all

doPlots = 0;  %Plot each road's source and recovered TSD signal. Keep disabled if you need to economize memory!
doLSPT = 0;   %Boolean -> Tell the code to optionally back-calculate the subgrade's loss of support under the joint.

%% update v2022-04-21
%for faster code, disable the " matrix close to signular" warnings that may
%pop up. They slow the code too much!
%Follow: https://www.mathworks.com/help/matlab/ref/lastwarn.html#responsive_offcanvas

load matlabWarningMessage
%this tiny variable has the warning message text and ID for the 'singular
%matrix warning" 
%turn off the warning
warnStruct = warning('off',warnID);
%   RESTORE IT AT THE END OF THE CODE!

%% addpath to dependencies (Denoising and back-calc engine)
addpath('back-calc\')
addpath('waveletDenoising\')%

%% EXPORT FILENAME
exportFilename = 'MnROAD_239_v20230425.xlsx';

%road id information
roadName = 'MnROAD LVR'
roadID = '239'

disp('Loading TSD Data')
inputFilename = 'T17202109270013_5cm.xlsx';
inputSheet    = 'TSD';

%%full LVR data range --> b2:b224011
%These calls are for SECTION 239 ONLY!

latLongData = xlsread(inputFilename,inputSheet,'b67603:c69003');       %%lat long data in decimal degrees
stationData = xlsread(inputFilename,inputSheet,'a67603:a69003');       %% station in m

loadData  = xlsread(inputFilename,inputSheet,'e67603:e69003'); %%dynamic load on the right-side half-axle [kg]!
loadData  = loadData.*9.81;  %parse to Newtons!

vxData = xlsread(inputFilename,inputSheet,'d67603:d69003');

vyData = xlsread(inputFilename,inputSheet,'g67603:p69003');
vyData = vyData(:,[7,6,5,4,3,2,1]);    %< NOTE THAT THIS COMES FROM THE TSD IN MM/SEC

TSDpoints = [0.130,0.215,0.300,0.450,0.600,0.900,1.500];
TSDpoints = TSDpoints';
numTSDSensors = length(TSDpoints);

%% THICKNESS DATA Received 2023-05-01/05. - this value is in inches!
thickness = 4.*ones(size(stationData));

%% 3) WAVELET DENOISING OF THE TSD DATA - 
disp('TSD denoising via wavelet decomposition')

latFrom = latLongData(:,1);
longFrom= latLongData(:,2);
clear latLongData

%% 3.0) Do the data denoising    
denoisedTSDvy  = zeros(size(vyData));

for zz = numTSDSensors:-1:1  %do this iteration back-wards on purpose so that I can plot the denoised SL110 with the "temporary names" from the denoising stage.
    refVY  = vyData(:,zz);
    fprintf('\t Denoising signal from TSD sensor at %g \n',TSDpoints(zz))        
    % 3.1: Data Denoising by wavelet decomposition.
    denoisedVy  = Haar_Denoise_LFDR(refVY,0.01);
    denoisedTSDvy(:,zz) = denoisedVy;
end  

clear refVY
clear denoisdedVy
%% 4) Do the plot of the recovered signal    
if doPlots        
    %%add an extra plot with all the denoised TSD signals       
    figure(300)
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);  %this should make the figure full-screen
    plotName = sprintf('Denoised TSD V_y data - %s, section %s',string(roadName), string(roadID));
    set(gcf,'Name',plotName);
    plot(stationData, denoisedTSDvy(:,1),'linewidth',1);
    set(gca,'FontSize',16)
    grid on
    titleString = sprintf('TSD signal for road %s, section %s. All denoised TSD sensors',string(roadName), string(roadID));
    title(titleString)
    xlabel ('station [km]')
    ylabel ('V_y [mm/s]')
    hold on
    for zz = 2:numTSDSensors
        plot(stationData, denoisedTSDvy(:,zz),'linewidth',1);        
    end        
%     legend('vy_{-450}','vy_{-300}','vy_{-200}','vy_{130}','vy_{210}','vy_{310}','vy_{450}','vy_{600}','vy_{900}', 'vy_{1510}')
    legend('vy_{130}','vy_{210}','vy_{310}','vy_{450}','vy_{600}','vy_{900}', 'vy_{1510}')
    hold off
    
end  %endif doPlots
drawnow

%% 4) Do the plot of the recovered signal    
if doPlots        
    %%add an extra plot with all the denoised TSD signals       
    figure(200)
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);  %this should make the figure full-screen
    plotName = sprintf('Source TSD V_y data- %s, section %s',string(roadName), string(roadID));
    set(gcf,'Name',plotName);
    plot(stationData, vyData(:,1),'linewidth',1);
    set(gca,'FontSize',16)
    grid on
    titleString = sprintf('TSD signal for road %s, section %s. All denoised TSD sensors',string(roadName), string(roadID));
    title(titleString)
    xlabel ('station [km]')
    ylabel ('V_y [mm/s]')
    hold on
    for zz = 2:numTSDSensors
        plot(stationData, vyData(:,zz),'linewidth',1);        
    end        
  %     legend('vy_{-450}','vy_{-300}','vy_{-200}','vy_{130}','vy_{210}','vy_{310}','vy_{450}','vy_{600}','vy_{900}', 'vy_{1510}')
    legend('vy_{130}','vy_{210}','vy_{310}','vy_{450}','vy_{600}','vy_{900}', 'vy_{1510}')
    hold off    
end  %endif doPlots
drawnow

%% 5) update v2022-02-28. Manual selector of weak joints to back-calculate
%     selectJoint = input('Do you want to back-calculate joints?. 1 if so, 0 to stop...');
selectJoint = 1;
while selectJoint
	%Automated joint search using FINDPEAKS on the SL 310
    %- using a minimum peak height of 75% quantile, and minimumdistance
    %between peaks of ~25-30 measurements seems to filter out all the
    %joint peaks (and ignore double peaks
    [~,index] = findpeaks(denoisedTSDvy(:,TSDpoints == 0.300),'MinPeakDistance',25,'MinPeakHeight',quantile(denoisedTSDvy(:,TSDpoints==0.300),0.75));
    jointPosition = stationData(index);   %% <<--- these are the station values at which the joints occur 
    jointLatLong = [latFrom(index), longFrom(index)];
		
	nnn = length(jointPosition);
    fprintf('%g joints were detected in this section \n',nnn)   

    %update v2023-0425. Solve the joint spacing to tell if the input
    %segment is a short pavement.
    jointSpacing = jointPosition(2:end) - jointPosition(1:end-1);
    jointSpacing = mean(jointSpacing);
    
    fprintf('The mean spacing between joints is %g \n',jointSpacing)   

    if jointSpacing.*1000 < 2.00 %meters, this tells about short slabs
        shortSlabs = 1;
    else
        shortSlabs = 0;
    end

	for i = 1:nnn		
        jointIndex = index(i);
        fprintf('processing joint at station [m] %g  \n',1000.*jointPosition(i));
        fprintf('\t solving continuous component \n')
		
        verborragia = 0; %put 1 if you want a verbose report of the back-calc on screen.
		
        subgradeType = 0;  % Boolean: 0 = Winkler foundation, 1 = Pasternak foundation.
        nu = 0.21;
        pressure = 115./145.04.*1e6; %TSD Wheel load -> 115 PSI to PA
            
        if ~shortSlabs
            %% back-calculation fro regular-sized slabs. Divide the problem in two as per TRB-2023 paper (Dissertation paper IV), with Aftermath04 update.
            %STAGE 1:
            %solve k,E,G ahead of the joint (mid-slab locations).
		    %DO MULTIPLE BACK-CALC OF THE CONTINUOUS COMPONENT FOR K, E, G, keep the mean value as representative one.
            stationAhead = stationData(max(jointIndex-45,1):max(jointIndex-35,1));  %Consider the deflection 2.25-1.75m ahead as 'mid-slab'
    
            localE = ones(size(stationAhead));
            localK = ones(size(stationAhead));
            localG = zeros(size(stationAhead));
            for jj = 1:length(stationAhead)
               
                localLoad = loadData(stationData==stationAhead(jj)).*1;  %leave in Newtons
                localThck = thickness(stationData==stationAhead(jj));  %% BY DEFAULT IT COMES IN INCHES!
                localThck = localThck.*2.54./100; %pass localThck to meters!            
                localvx   = vxData(stationData==stationAhead(jj));
                localTSD_short = denoisedTSDvy(stationData==stationAhead(jj),:);
                
                %update v04-21: Call the back-calculation front-end based on deflection velocity!
                %% careful here! localTSD_cont is in mm/sec. Must pass it to the back-calc engine in m/sec!!!!
                %% also, localThick is in inches, must pass to meters!
                %% stability update v2023-03-30.
                %Move the NaN check to here because beforehand, I don't know
                %which sensor may drop a NaN
                %Do the check now, and if I find nans, add/remove sensors
                %accordingly.
                notNANmeas = ~isnan(localTSD_short);
                localTSD_short = localTSD_short(notNANmeas);
                localTSDpoints = TSDpoints(notNANmeas);
    
                [localE(jj),localK(jj),localG(jj),~,~,~,~] =  backCalc_continuous_0420(localTSDpoints,localTSD_short./1e3,localvx,localLoad,pressure,localThck,nu,subgradeType,verborragia);      
            end
           
            %get the final k, E, G as definitivevalues
            localE = mean(localE);
            localK = mean(localK);
            localG = mean(localG);
            
            %STAGE 2 - updated v2023-01-25
            %DO THE BACK-CALC BASED ON DEFLECTION SPEED. 
            %NO NEED TO CHANGE SIGNS (TSD'S convention on vy is the same as
            %for my w(z,t), a pavement that goes down has positive vy
            fprintf('\t solving joint approach \n')
            stationRanges = max(jointIndex-30,1):1:max(jointIndex-4,1);  %% <-- all these positions are the ones I'm back-calculating, start 1.5m ahead, stop 20cm ahead of the joint.
		    nk = length(stationRanges);             
            
		    jointLocation = zeros(nk,1);
		    LTE = zeros(nk,1);
		    LossSupport = zeros(nk,1);
            SSEfinal = zeros(nk,1);
		    localTSD_short = zeros(nk,length(TSDpoints));
		    
            for j = 1:nk
                localTSD_short(j,:) = denoisedTSDvy(stationRanges(j),:); 
                localThck = thickness(stationRanges(j));
                localThck = localThck.*2.54./100; %pass localThck to meters!
                localLoad = loadData(stationRanges(j)).*1;  %leave in Newtons
                % 
                jointLocation(j) = -1;
                LTE(j) = -1;
                LossSupport(j) = 1;        
                fprintf('\t solving joint profile at station [m] %g \n',1000.*stationData(stationRanges(j)))
                %update V2023-01-25: Get an estimate of the distance between TSD wheel and joint, use that for the back-calc
                %for LTE. 
                jointLocation(j) = jointPosition(i) - stationData(stationRanges(j));  %The distance between the current measurement and the joint's location should be a rough estimate of where the joint is at [+/- 10-15cm]
                %--> safety check: if jointLocation(j) is equal to any of
                %the TSD sensors locations (where vy goes to infinity),
                % shift it by 2cm (add error to prevent divergence)
                estimateCSafetyCheck = find(TSDpoints == jointLocation(j));
                if ~isempty(estimateCSafetyCheck)
                    %the joint location is falling over a TSD sensor,
                    %assume it elsewhere (shift +1cm) to prevent divergence
                    jointLocation(j) = jointLocation(j)+0.01;
                end        
                %% solving the joint
                %% IMPORTANT: VYcomes in mm/s, and vx in m/s. PASS BOTH IN M/sec	        
                % update v2023-01-25 -> Remove correction for estimateC,
                % solve for LTE and Loss of support only!         
				
                %Update V 2022-03-23 -> Use back-calc based on deflection velocity! 
                % %May need vy and local vx record too!
                localvx = vxData(stationRanges(j));
                        
                %stability update v2023-03-30 -> add the nan-check here!
                notNANmeas = ~isnan(localTSD_short(j,:));
                good_localTSD_pulse = localTSD_short(j,notNANmeas);
                localTSDpoints = TSDpoints(notNANmeas);	        
                [LTE(j),LossSupport(j),SSEfinal(j),~,~,~] =  backCalc_joint_0420(localTSDpoints,good_localTSD_pulse./1e3,localvx,jointLocation(j),localThck,localK,localG,localE,nu,localLoad,pressure,doLSPT,verborragia);
                fprintf('\t joint"s LTE is %g \n',LTE(j))
            end

        else
            %update v2023-04-25 -> Special back-calc solver for short
            %slabs. Solve k,E,G, LTE (disregard LSPT)

            stationRanges = max(jointIndex-30,1):1:max(jointIndex-4,1);  %% <start back-calc 1.5m ahead, stop 20cm ahead of the joint.
		    nk = length(stationRanges);             
            
		    jointLocation = zeros(nk,1);
            localK = zeros(nk,1);
            localE = zeros(nk,1);
            localG = zeros(nk,1);
		    LTE = zeros(nk,1);
%             LossSupport = zeros(nk,1); 
            
		    SSEfinal = zeros(nk,1);
		    
            for j = 1:nk
                localThck = thickness(stationRanges(j));
                localThck = localThck.*2.54./100; %pass localThck to meters!
                localLoad = loadData(stationRanges(j)).*1;  %leave in Newtons
                localvx = vxData(stationRanges(j));
                
                fprintf('\t solving joint profile at station [m] %g \n',1000.*stationData(stationRanges(j)))
                
                %The distance between the current measurement and the joint's location should be a rough estimate of where the joint is at [+/- 10-15cm]
                jointLocation(j) = jointPosition(i) - stationData(stationRanges(j)); 
                
                %get the local TSD vy basin and discard any nans
                localTSD_short = denoisedTSDvy(stationRanges(j),:); 
                notNANmeas = ~isnan(localTSD_short);
                localTSD_short = localTSD_short(notNANmeas);
                localTSDpoints = TSDpoints(notNANmeas);
    
                [localE(j),localK(j),localG(j),LTE(j),SSEfinal(j)] = backCalc_shortSlab(localTSDpoints,localTSD_short./1e3,localvx,jointLocation(j),localLoad,pressure,localThck,nu,subgradeType,verborragia);      
                fprintf('\t joint"s LTE is %g \n',LTE(j))
            end
            %for the fun of it, estimate loss of support over k-value by
            %comparing the kValue at the 1.50m measurement to the remaining
            %ones
            LossSupport = localK./(localK(1));
            
        end
		
		%% update v 2022-03-04. Do the export of the results now. I don't want to deal with the structure (that can be messy to export in few steps)
        %update v2023-03-31 -> Since there're too many joints in this case
        %study, store them in tiny matlab files.
        
		exportMini = sprintf('joint_station_%g.mat',jointPosition(i));
		
        %%update v2023-03-30 -> Export summary results to Excel (saves time)
        % exportFilename stated at the beginning!
		exportRow = 4 + i;
        exportRange = sprintf('b%g:j%g',exportRow,exportRow);
        %update v2023-04-25 -> Force a mean to the kEG back-calcs here
        %because the short-slab procedure outputs vector data and not a
        %single value!
		exportVariable = [jointPosition(i),jointLatLong(i,:),mean(localK),mean(localE),mean(localG),mean(LTE),mean(LossSupport),mean(SSEfinal)];
		export = xlswrite(exportFilename,exportVariable,'summary',exportRange);
		toc
	end  %end for i = 1 : nnn
	%% before closing the while, check if doing one more joint.
	%selectJoint = input('Joint Completed. Solve another one?. 1 if so, 0 to stop...');
	selectJoint = 0;
end
save MnRoad_239_v20230425.mat

disp('....all completed')
toc
%% RESTORE THE SINGULAR MATRIX WARNING
warning(warnStruct);
