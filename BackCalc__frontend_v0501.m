%% Concrete back-calculation from TSD: FRONT-END SCRIPT FOR LTE BACK-CALCULATION
%% TSD 5-cm dataset from MnROAD section #239

%V6.1 - 2024-01-26
% --Updated output file: add header row to the output table describing each column.
%V6.0 - 2023-04-25 
% -- Updated front=end to solve the case of pavements with short slabs
% (back-calculation that simultaneously solves the k,E,G and the LTE of the % joint)
% -- added numerical solver for LTE: removed correction for 'c'
% -- Experimental: Added estimate of subgrade loss of support (LSPT). 
% User must turn variable doLSPT = 1 to do the calculation.

%% STAGE 0 DATA LOAD-----
%LOAD HERE THE TSD RECORDS

tic
clc

restoredefaultpath
clear variables
close all

doPlots = 0;  		%Boolean: 1 = Plot each road's input and denoised TSD signal. Keep disabled if you need to economize memory!
doLSPT = 0;   		%Boolean: [EXPERIMENTAL] 1 = Tell the code to optionally back-calculate the subgrade's loss of support under the joint.
verborragia = 0;	%Boolean: 1 = if you want a verbose report of the back-calc progress on screen.
subgradeType = 0;  	%Boolean: 0 = Winkler foundation, 1 = Pasternak foundation.
doBackCalc = 1;		%Boolean: 1 = do the back-calculation process. 0 = limit to import and denoise the TSD vy data.
    
%% update v2022-04-21
%For faster code, disable the "matrix close to singular" warnings that may pop up. They slow the code too much!
%Follow: https://www.mathworks.com/help/matlab/ref/lastwarn.html#responsive_offcanvas

load matlabWarningMessage;  %<---this tiny variable has the warning message text and ID for the 'singular matrix warning'
%Turn off the warning
warnStruct = warning('off',warnID);
% <<RESTORE IT AT THE END OF THE CODE!>>

%% addpath to dependencies (Denoising and back-calc engine)
addpath('back-calc\')
addpath('waveletDenoising\')

%% EXPORT FILENAME. User: change this output file name as you see fit. 
exportFilename = 'MnROAD_239_v20230425.xlsx';

%% ROAD ID Information. User: change this output file name as you see fit. 
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
vyData = xlsread(inputFilename,inputSheet,'g67603:p69003');   %<Units: MM/SEC

%% tiny re-sorting of vy data so that they match the location of the Doppler sensors as per the "TSDpoints" variable below 
vyData = vyData(:,[7,6,5,4,3,2,1]);   
%vyData = vyData(:,[10,9,8,7,6,5,4,3,2,1]);    

TSDpoints = [0.130,0.215,0.300,0.450,0.600,0.900,1.500];
%TSDpoints = [-0.45,-0.300,-0.200,0.130,0.215,0.300,0.450,0.600,0.900,1.500];
TSDpoints = TSDpoints';
numTSDSensors = length(TSDpoints);

%% THICKNESS DATA. Loaded value is in inches. Convert to meters.
thickness = 4.*ones(size(stationData));
thickness = thickness .*2.54./100;

latFrom = latLongData(:,1);
longFrom= latLongData(:,2);
clear latLongData

%% 3) WAVELET DENOISING OF THE TSD DATA 
disp('TSD denoising via Haar wavelet decomposition with local adaptive FDR [Katicha et al., 2024]')

%% 3.0) Do the data denoising    
denoisedTSDvy  = zeros(size(vyData));

for zz = 1:numTSDSensors
    refVY  = vyData(:,zz);
    fprintf('\t Denoising signal from TSD sensor at %g \n',TSDpoints(zz))        
    
    denoisedVy  = Haar_Denoise_LFDR(refVY,0.01);
    denoisedTSDvy(:,zz) = denoisedVy;
end  

clear refVY
clear denoisdedVy

%% 4) Do the plot of the recovered signal    
if doPlots        
    %%Plot all the denoised TSD signals       
    figure(300)
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);  %this should make the figure full-screen
    plotName = sprintf('Denoised TSD V_y data - %s, section %s',string(roadName), string(roadID));
    set(gcf,'Name',plotName);
    plot(stationData, denoisedTSDvy(:,1),'linewidth',1);
    set(gca,'FontSize',16)
    grid on
    titleString = sprintf('TSD signal for road %s, section %s. All denoised TSD sensors',string(roadName), string(roadID));
    title(titleString)
    xlabel ('station [m]')
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

%% 4) Do the plot of the original noisy signal    
if doPlots        
    figure(200)
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);  %this should make the figure full-screen
    plotName = sprintf('Source TSD V_y data- %s, section %s',string(roadName), string(roadID));
    set(gcf,'Name',plotName);
    plot(stationData, vyData(:,1),'linewidth',1);
    set(gca,'FontSize',16)
    grid on
    titleString = sprintf('TSD signal for road %s, section %s. All noisy TSD sensors',string(roadName), string(roadID));
    title(titleString)
    xlabel ('station [m]')
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

%% DO THE BACK-CALCULATION!
%% 5) update v2022-02-28. Manual selector of weak joints to back-calculate
while doBackCalc
    %Automated joint search using FINDPEAKS on the SL 310 using a minimum peak height of 75% quantile,
    % and minimum distance between peaks of ~25-30 measurements seems to filter out all the
    %joint peaks (and ignore double peaks
    
    [~,index] = findpeaks(denoisedTSDvy(:,TSDpoints == 0.300),'MinPeakDistance',25,'MinPeakHeight',quantile(denoisedTSDvy(:,TSDpoints==0.300),0.75));
    jointPosition = stationData(index);   %% <<--- these are the station values at which the joints occur 
    jointLatLong = [latFrom(index), longFrom(index)];
    
    nnn = length(jointPosition);
    fprintf('%g joints were detected in this section \n',nnn)   

    %update v2023-0425. Solve the mean spacing between joints to tell if the input segment is a short pavement.
    %IMPORTANT -- ENSURE THAT JOINT-POSITION (AND STATION) ARE IN METERS WHEN YOU GET TO THIS POINT!
    
    jointSpacing = jointPosition(2:end) - jointPosition(1:end-1);
    jointSpacing = mean(jointSpacing);    
    fprintf('The mean spacing between joints is %g \n',jointSpacing)   

    if jointSpacing < 2.00 %meters, this tells about short slabs
        shortSlabs = 1;
    else
        shortSlabs = 0;
    end
    nu = 0.21;
    pressure = 115./145.04.*1e6; %TSD Wheel load -> 115 PSI to PA
    
    for i = 1:nnn		
        jointIndex = index(i);
        fprintf('processing joint at station [m] %g  \n',jointPosition(i));
        fprintf('\t solving continuous component \n')
            
        if ~shortSlabs
            %% back-calculation foR regular-sized slabs [slab length >2.00m]. Divide the problem in two as per TRB-2023 paper (Dissertation paper IV), with Aftermath04 update.
            %STAGE 1:
            %solve k,E,G ahead of the joint (mid-slab locations).
	    %DO MULTIPLE BACK-CALC OF THE CONTINUOUS COMPONENT FOR K, E, G, keep the mean value as representative one.
            stationAhead = stationData(max(jointIndex-45,1):max(jointIndex-35,1));  %Consider the deflection 2.25-1.75m ahead as 'mid-slab'
    
            localE = ones(size(stationAhead));
            localK = ones(size(stationAhead));
            localG = zeros(size(stationAhead));
            for jj = 1:length(stationAhead)
               
                localLoad = loadData(stationData==stationAhead(jj)).*1; %Units: Newtons
                localThck = thickness(stationData==stationAhead(jj));   %Update v2024-01-26. This is now in meters!                        
                localvx   = vxData(stationData==stationAhead(jj));
                localTSD_short = denoisedTSDvy(stationData==stationAhead(jj),:);
                
                %update v04-21: Call the back-calculation front-end based on deflection velocity!
                %% IMPORTANT! localTSD_cont is in mm/sec. Must pass it to the back-calc engine in m/sec!!!!
                
		%% stability update v2023-03-30.
                %Move the NaN check to here because beforehand, I don't know which sensor may drop a NaN
                %Do the check now, and if I find nans, add/remove sensors accordingly.
		
                notNANmeas = ~isnan(localTSD_short);
                localTSD_short = localTSD_short(notNANmeas);
                localTSDpoints = TSDpoints(notNANmeas);
    
                [localE(jj),localK(jj),localG(jj),~,~,~,~] =  backCalc_continuous_0420(localTSDpoints,localTSD_short./1e3,localvx,localLoad,pressure,localThck,nu,subgradeType,verborragia);      
            end
           
            %get the final k, E, G as definitive values
            localE = mean(localE);
            localK = mean(localK);
            localG = mean(localG);
            
            %STAGE 2 - updated v2023-01-25
            %DO THE BACK-CALC BASED ON DEFLECTION SPEED. 
            
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
                localThck = thickness(stationRanges(j));    %Update 2024-01-26: This is in meters now!               
                localLoad = loadData(stationRanges(j)).*1;  %leave in Newtons
		 
                jointLocation(j) = -1;
                LTE(j) = -1;
                LossSupport(j) = 1;        
                fprintf('\t solving joint profile at station [m] %g \n',stationData(stationRanges(j)))
                %update V2023-01-25: Get an estimate of the distance between TSD wheel and joint [Variable jointLocation]
		%Use that for the back-calc for LTE. 
                jointLocation(j) = jointPosition(i) - stationData(stationRanges(j));  %The distance between the current measurement and the joint's location should be a rough estimate of where the joint is at [+/- 10-15cm]
               
		%--> safety check: if jointLocation(j) is equal to any of the TSD sensors locations (where vy goes to infinity),
                % shift it by 2cm (add error to prevent divergence)
                estimateCSafetyCheck = find(TSDpoints == jointLocation(j));
                if ~isempty(estimateCSafetyCheck)
                    %the joint location is falling over a TSD sensor,
                    %assume it elsewhere (shift +1cm) to prevent divergence
                    jointLocation(j) = jointLocation(j)+0.01;
                end  
		
                %% Solving the joint's LTE
                %% IMPORTANT: VY comes in mm/s, and vx in m/s. PASS BOTH IN M/sec	        
                % update v2023-01-25 -> Remove correction for estimateC,
                % solve for LTE and Loss of support only!         
				
                %Update V 2022-03-23 -> Use back-calc based on deflection velocity! 
                %Needs vy and local vx record too!
                localvx = vxData(stationRanges(j));
                        
                %stability update v2023-03-30 -> add the nan-check here!
                notNANmeas = ~isnan(localTSD_short(j,:));
                good_localTSD_pulse = localTSD_short(j,notNANmeas);
                localTSDpoints = TSDpoints(notNANmeas);	        
                [LTE(j),LossSupport(j),SSEfinal(j),~,~,~] =  backCalc_joint_0420(localTSDpoints,good_localTSD_pulse./1e3,localvx,jointLocation(j),localThck,localK,localG,localE,nu,localLoad,pressure,doLSPT,verborragia);
                fprintf('\t joint""s LTE is %g \n',LTE(j))
            end
        else
            %update v2023-04-25: Special back-calc solver for short slabs [slab length < 2.00m]. 
	    %Solve k, E, G, LTE (disregard LSPT)

            stationRanges = max(jointIndex-30,1):1:max(jointIndex-4,1);  %% <start back-calc 1.5m ahead, stop 20cm ahead of the joint.
	    nk = length(stationRanges);  
            
	    jointLocation = zeros(nk,1);
            localK = zeros(nk,1);
            localE = zeros(nk,1);
            localG = zeros(nk,1);
	    LTE = zeros(nk,1);
            SSEfinal = zeros(nk,1);
	    %LossSupport = zeros(nk,1);   %% no need to get LSPT in this solver, as k is corrected as I approach the joint.
            	    
            for j = 1:nk
                localThck = thickness(stationRanges(j)); 	%Update v2024-01-26: localThck now in meters!
                localLoad = loadData(stationRanges(j)).*1;  	%Units: Newtons
                localvx = vxData(stationRanges(j));
                
                fprintf('\t solving joint profile at station [m] %g \n',stationData(stationRanges(j)))
                
                %The distance between the current measurement and the joint's location should be a rough estimate of where the joint is at [+/- 10-15cm]
                jointLocation(j) = jointPosition(i) - stationData(stationRanges(j)); 
                
                %get the local TSD vy basin and discard any nans
                localTSD_short = denoisedTSDvy(stationRanges(j),:); 
                notNANmeas = ~isnan(localTSD_short);
                localTSD_short = localTSD_short(notNANmeas);
                localTSDpoints = TSDpoints(notNANmeas);
    
                [localE(j),localK(j),localG(j),LTE(j),SSEfinal(j)] = backCalc_shortSlab(localTSDpoints,localTSD_short./1e3,localvx,jointLocation(j),localLoad,pressure,localThck,nu,subgradeType,verborragia);      
                fprintf('\t joint""s LTE is %g \n',LTE(j))
            end
            %for the fun of it, estimate the loss of support over k-value by
            %comparing the k - Value at the 1.50m measurement to the remaining ones
            LossSupport = localK./(localK(1));            
        end
		
	%%Update v2022-03-04. Do the export of the results now.
        %%update v2023-03-30. Export summary results to Excel (saves time)
	%%Update v2024-01-26. Add header row to Excel Export
 	
  	if i == 1
            header = {'station [m]','lat','long','k [N/m3]','E [N/m2]','G [N/m]','LTE','LSPT','SSE'};
            xlswrite(exportFilename,header,exportSheet,'b3:j3');
        end	
 
        % exportFilename stated at the beginning!
	exportRow = 4 + i;
        exportRange = sprintf('b%g:j%g',exportRow,exportRow);
	
	%update v2023-04-25: Force a mean to the k,E,G back-calcs here
        %because the short-slab procedure outputs vector data and not a single value!
		exportVariable = [jointPosition(i),jointLatLong(i,:),mean(localK),mean(localE),mean(localG),mean(LTE),mean(LossSupport),mean(SSEfinal)];
		export = xlswrite(exportFilename,exportVariable,'summary',exportRange);
		toc
	end  %end for i = 1 : nnn
	%% Finished processing joints, exit the 'while' loop.
	doBackCalc = 0;
end

%optional - save progress in a big Mat file. USER: Change the name of the output file as you see fit.
%save MnRoad_239_v20230425.mat

disp('....all completed')
toc

%% Restore the "singular matrix" warning.
warning(warnStruct);
