%% DISSERTATION PAPER 3 - FRONT-END SCRIPT FOR TSD BACK-CALCULATION
%%- TSD VALIDATION RUN 5cm from MnROAD Loop section 124
%Candidate Release V2022-05-01

%% STAGE 0 DATA LOAD
%LOAD HERE THE TSD RECORDS
tic
clc
restoredefaultpath
clear variables
close all
doPlots = 1;  %Plot each road's source and recovered TSD signal. Keep disabled if you need to economize memory!

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
% ....

%% addpath to the V2.0-compatible dependencies! (Denoising and back-calc engine)
addpath('./waveletDenoising')
addpath('./back-calc')

% 
disp('Loading TSD Data')
inputFilename = 'T17202109270013_5cm.xlsx';
inputSheet    = 'TSD';

%% input block for LVR sections 124-524
% latLongData = xlsread(inputFilename,inputSheet,'b29823:c32723');       %%lat long data in decimal degrees
%section 124
stationData = xlsread(inputFilename,inputSheet,'a29823:a32723');       %% station in m
roadID = 124.*ones(size(stationData));
roadName = roadID ; 
loadData  = xlsread(inputFilename,inputSheet,'e29823:e32723'); %%dynamic load on the right-side half-axle [kg]!
loadData  = loadData.*9.81;  %parse to Newtons!
loadInLBS = 0; %<<--put 1 if the load is in LBS. Otherwise, assume NEWTONS.
%WE KNOW THE CONCRETE SLABS ARE 6-INCHES ALL THROUGHOUT THE SECTION. BUILT 2017.
%NO MATERIAL TEST INFO ON THE MNROAD DATABASE. LTE TESTED BY FWD IN 2019
lyrThick = 6.*ones(size(stationData));
%The export will generate a sheet with the processed joint info once it completes a loop in the 'while' statement.
exportFilename = 'MNRoad_Loop524_5cm_v0420.xlsx';
% update v2022-04-20 -> Use deflection velocity data only!

%columns are::  Slope 1.500[µm/m]	 Slope 0.900[µm/m]	 Slope 0.600[µm/m]	
%               Slope 0.450[µm/m]	 Slope 0.300[µm/m]	 Slope 0.215[µm/m]	 Slope 0.130[µm/m]	
%               Slope -0.200[µm/m]	 Slope -0.300[µm/m] Slope -0.450[µm/m]
vyData = xlsread(inputFilename,inputSheet,'g29823:p32723');
vxData = xlsread(inputFilename,inputSheet,'d29823:d32723');

vyData = vyData(:,[10,9,8,7,6,5,4,3,2,1]);    %< NOTE THAT THIS COMES FROM THE TSD IN MM/SEC
TSDpoints = [-0.45,-0.3,-0.2,0.130,0.215,0.300,0.450,0.600,0.900,1.500];
TSDpoints = TSDpoints';
numTSDSensors = length(TSDpoints);

stationFrom = stationData(:,1);  %%get station information (from and to) for each measurement.
% latFrom = latLongData(:,1);
% longFrom= latLongData(:,2);
% latTo   = latFrom;
% longTo  = longFrom;


%% input block for section 138-238
%section 238
% stationData = xlsread(inputFilename,inputSheet,'a63003:a65603');       %% station in m
% roadID = 238.*ones(size(stationData));
% roadName = roadID ; 
% 
% loadData  = xlsread(inputFilename,inputSheet,'e63003:e65603'); %%dynamic load on the right-side half-axle [kg]!
% loadData  = loadData.*9.81;  %parse to Newtons!
% loadInLBS = 0; %<<--put 1 if the load is in LBS. Otherwise, assume NEWTONS.
% 
% %WE KNOW THE CONCRETE SLABS ARE 6-INCHES ALL THROUGHOUT THE SECTION. BUILT 2017.
% %NO MATERIAL TEST INFO ON THE MNROAD DATABASE. LTE TESTED BY FWD IN 2019
% lyrThick = 8.*ones(size(stationData));
% 
% exportFilename = 'MNRoad_Loop238_5cm_v0420.xlsx';
% 
% %columns are::  Slope 1.500[µm/m]	 Slope 0.900[µm/m]	 Slope 0.600[µm/m]	
% %               Slope 0.450[µm/m]	 Slope 0.300[µm/m]	 Slope 0.215[µm/m]	 Slope 0.130[µm/m]	
% %               Slope -0.200[µm/m]	 Slope -0.300[µm/m] Slope -0.450[µm/m]
% vyData = xlsread(inputFilename,inputSheet,'g63003:p65603');
% vxData = xlsread(inputFilename,inputSheet,'d63003:d65603');
% 
% vyData = vyData(:,[10,9,8,7,6,5,4,3,2,1]);    %< NOTE THAT THIS COMES FROM THE TSD IN MM/SEC
% TSDpoints = [-0.45,-0.3,-0.2,0.130,0.215,0.300,0.450,0.600,0.900,1.500];
% TSDpoints = TSDpoints';
% numTSDSensors = length(TSDpoints);
% 
% stationFrom = stationData(:,1);  %%get station information (from and to) for each measurement.
% % latFrom = latLongData(:,1);
% % longFrom= latLongData(:,2);
% % latTo   = latFrom;
% % longTo  = longFrom;

%% input block for section 239
%section 239
% stationData = xlsread(inputFilename,inputSheet,'a67703:a68803');       %% station in m
% roadID = 239.*ones(size(stationData));
% roadName = roadID ; 
% 
% loadData  = xlsread(inputFilename,inputSheet,'e67703:e68803'); %%dynamic load on the right-side half-axle [kg]!
% loadData  = loadData.*9.81;  %parse to Newtons!
% loadInLBS = 0; %<<--put 1 if the load is in LBS. Otherwise, assume NEWTONS.
% 
% %WE KNOW THE CONCRETE SLABS ARE 6-INCHES ALL THROUGHOUT THE SECTION. BUILT
% %2017. NO MATERIAL TEST INFO ON THE MNROAD DATABASE. LTE TESTED BY FWD IN
% %2019
% lyrThick = 4.*ones(size(stationData));
% 
% %columns are::  Slope 1.500[µm/m]	 Slope 0.900[µm/m]	 Slope 0.600[µm/m]	
% %               Slope 0.450[µm/m]	 Slope 0.300[µm/m]	 Slope 0.215[µm/m]	 Slope 0.130[µm/m]	
% %               Slope -0.200[µm/m]	 Slope -0.300[µm/m] Slope -0.450[µm/m]
% vyData = xlsread(inputFilename,inputSheet,'g67703:p68803');
% vxData = xlsread(inputFilename,inputSheet,'d67703:d68803');
% 
% %vyData = vyData(:,[10,9,8,7,6,5,4,3,2,1]);    %< NOTE THAT THIS COMES FROM THE TSD IN MM/SEC
% %TSDpoints = [-0.45,-0.3,-0.2,0.130,0.215,0.300,0.450,0.600,0.900,1.500];
% 
% %Short-slab case study: Use without trailing sensors!
% vyData = vyData(:,[7,6,5,4,3,2,1]);    %< NOTE THAT THIS COMES FROM THE TSD IN MM/SEC
% TSDpoints = [0.130,0.215,0.300,0.450,0.600,0.900,1.500];
% 
% TSDpoints = TSDpoints';
% numTSDSensors = length(TSDpoints);
% 
% stationFrom = stationData(:,1);  %%get station information (from and to) for each measurement.
% % latFrom = latLongData(:,1);
% % longFrom= latLongData(:,2);
% % latTo   = latFrom;
% % longTo  = longFrom;
% 

%% Get all unique roads in the TSD dataset - 
% UPDATE V2019-10-30: MUST DISTINGUISH BY roadName+blockID [roadID and roadNAME by themselves would collect many segments together (and overlap signals)]
roadIDList = unique(roadID);    %use this vector to locate the TSD points for each road.
numRoads = length(roadIDList);
fprintf('\t %g streets (+ blocks) were recognized \n',numRoads)

%% 3) WAVELET DENOISING OF THE TSD DATA - 
%%update v 2022-03-17 St Patrick -> Do Haar Wavelet decomposition only!
disp('TSD denoising via wavelet decomposition')

for k = 1:numRoads
    whichRoad = roadIDList(k);
%     tsdIndices = find(strcmp(roadID,whichRoad));  %this tweak allows to use find to compare cell arrays. Source: https://www.mathworks.com/matlabcentral/answers/84242-find-in-a-cell-array
    tsdIndices = find(roadID ==whichRoad);  
    %get the name of the road/lane under analysis
    thatRoadName = roadName(tsdIndices);
    thatRoadName = thatRoadName(1);            
    
    %get the TSD data for that road
    fprintf(' \t now processing: %s, segment %s \n',string(thatRoadName),string(whichRoad));
    
    vy = vyData(tsdIndices,:);
    vx = vxData(tsdIndices,:);
    loadRightWheel = loadData(tsdIndices);
    
    station = stationFrom(tsdIndices);
%     latStart = latFrom(tsdIndices);
%     latEnd   = latTo(tsdIndices);
%     longStart= longFrom(tsdIndices);
%     longEnd  = longTo(tsdIndices);
    thickness= lyrThick(tsdIndices);
    
    %Remove Nans that may cause trouble. Use this workaround: https://www.mathworks.com/matlabcentral/answers/164316-select-everything-not-returned-by-index    
    theseAreNans = find(isnan(vy(:,1))); %this returns a vector with the positions of NaN's inside the defl.Slope 110.
    NansIndex = isnan(vy(:,1));  %this logical array (different from the one above) helps then get the non-NaN values easily. [it's a 0/1 vector as long as refTSD in which the 1's indicate the NaNs]
    
    vy = vy(~NansIndex,:);
    thickness = thickness(~NansIndex,:);
    loadRightWheel = loadRightWheel(~NansIndex,:);
    station = station(~NansIndex,:);
    
    %% 3.0) Do the data denoising    
    denoisedTSDvy     = zeros(size(vy));
    
    for zz = numTSDSensors:-1:1  %do this iteration back-wards on purpose so that I can plot the denoised SL110 with the "temporary names" from the denoising stage.
        refVY  = vy(:,zz);
        fprintf('\t Denoising signal from TSD sensor at %g \n',TSDpoints(zz))        
        % 3.1: Data Denoising by wavelet decomposition.
        denoisedVy  = Haar_Denoise_LFDR(refVY,0.01);
        denoisedTSDvy(:,zz) = denoisedVy;
    end  
    %% 4) Do the plot of the recovered signal    
    if doPlots        
        figure(300+k)
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);  %this should make the figure full-screen
        plotName = sprintf('Denoised TSD All Channels - %s, section %s',string(thatRoadName), string(whichRoad));
        set(gcf,'Name',plotName);
        plot(station, denoisedTSDvy(:,1),'linewidth',1);
        set(gca,'FontSize',18)
        grid on
        titleString = sprintf('TSD signal for road %s, section %s. All denoised TSD sensors',string(thatRoadName), string(whichRoad));
        title(titleString)
        xlabel ('station [m]')
        ylabel ('vertical defl velocity [mm/s]')
        hold on
        for zz = 2:numTSDSensors
            plot(station, denoisedTSDvy(:,zz),'linewidth',1);        
        end        
        legend('vy_{-450}','vy_{-300}','vy_{-200}','vy_{130}','vy_{210}','vy_{310}','vy_{450}','vy_{600}','vy_{900}', 'vy_{1510}')
        hold off
        
    end  %endif doPlots
    drawnow
    
    %% 5) update v2022-02-28. Manual selector of weak joints to back-calculate
%     selectJoint = input('Do you want to back-calculate joints?. 1 if so, 0 to stop...');
    selectJoint = 1;
    jointsInRoadK = 0;
    while selectJoint
% 		jointPosition = input('give the joint"s station (where the TSD bumps) [m]...');
% 		jointPosition = [1551, 1555.8, 1560.5];
%       jointPosition = sort(jointPosition);
		        
		%experimental automated joint search using FINDPEAKS on the SL 310
        %- using a minimum peak height of 75% quantile, and minimumdistance
        %between peaks of ~25-30 measurements seems to filter out all the
        %joint peaks (and ignore double peaks
        [~,index] = findpeaks(denoisedTSDvy(:,TSDpoints == 0.300),'MinPeakDistance',25,'MinPeakHeight',quantile(denoisedTSDvy(:,TSDpoints==0.300),0.75));
        jointPosition = station(index);
        nnn = length(jointPosition);        
        sprintf('%g joints were detected in this section \n',nnn)        
               
		for i = 1:nnn
			jointIndex = find(station==jointPosition(i));
			%         jointLatLong = [latStart(jointPosition(i));longStart(jointPosition(i))];
			jointsInRoadK = jointsInRoadK+1;			
			fprintf('\t solving continuous component \n')
			subgradeType = 0;
            nu = 0.21;
            
            %% Solve k,E,G ahead of the joint.
			%UPDATE V 2022-03-23 -> DO MULTIPLE BACK-CALC OF THE CONTINUOUS
            %COMPONENT FOR K, E, G. keep the mean value as representative one.
            stationAhead = station(max(jointIndex-45,1):max(jointIndex-35,1));  %consider the deflection 2.25-1.75m ahead as 'mid-slab'
            localE = zeros(size(stationAhead));
            localK = zeros(size(stationAhead));
            localG = zeros(size(stationAhead));
            for jj = 1:length(stationAhead)
                if loadInLBS
                    localLoad = loadRightWheel(station==stationAhead(jj)).*4.44822162;%< pass axleLoad from LBS to Newton
                    localLoad = localLoad.*4.44822162;%< pass axleLoad from LBS to Newton    
                else
                    localLoad = loadRightWheel(station==stationAhead(jj)).*1;  %leave in Newtons
                    localLoad = localLoad.*1;%< leave in  Newton    
                end
                pressure = 115./145.04.*1e6; %TSD Wheel load -> 115 PSI to PA
                localThck = thickness(station==stationAhead(jj));
                localThck = localThck.*2.54./100;
                
                localvx = vx(station==stationAhead(jj));
                localTSD_cont = denoisedTSDvy(station==stationAhead(jj),:);
                
                %update v04-21: Call the back-calculation front-end based on deflection velocity!
                %% careful here! localTSD_cont is in mm/sec. Must pass it to the back-calc engine in m/sec!!!!
                %% also, localThick is in inches, must pass to meters!
                verborragia = 0;                
                [localE(jj),localK(jj),localG(jj),~,~,~,~] =  backCalc_continuous_0420(TSDpoints,localTSD_cont./1e3,localvx,localLoad,pressure,localThck,nu,subgradeType,verborragia);      
            end
            
            %get the final k, E, G as definitivevalues
            localE = mean(localE);
            localK = mean(localK);
            localG = mean(localG);
            
            %% Now proceed to the joint.
            %UPDATE V 2022-03-23: DO THE BACK-CALC BASED ON DEFLECTION SPEED. 
            %NO NEED TO CHANGE SIGNS (TSD'S convention on vy is the same as
            %for my w(z,t), a pavement that goes down has positive vy
            fprintf('\t solving joint approach \n')
            stationRanges = max(jointIndex-30,1):1:max(jointIndex-4,1);  %% <-- all these positions are the ones I'm back-calculating, stop 20cm ahead of the joint.
			nk = length(stationRanges);       
            
			jointLocation = zeros(nk,1);
			LTE = zeros(nk,1);
			LTE2 = zeros(nk,1);
            SSEfinal = zeros(nk,1);
			localTSD_pulse = zeros(nk,length(TSDpoints));
			
			for j = 1:nk
				localTSD_pulse(j,:) = denoisedTSDvy(stationRanges(j),:); 
				localThck = thickness(stationRanges(j));
				localThck = localThck.*2.54./100; %pass localThck to meters!
                
				if loadInLBS
					localLoad = loadRightWheel(stationRanges(j)).*4.44822162;%< pass axleLoad from LBS to Newton
					localLoad = localLoad.*4.44822162;%< pass axleLoad from LBS to Newton  
				else
					localLoad = loadRightWheel(stationRanges(j)).*1;  %leave in Newtons
					localLoad = localLoad.*1;%< leave in  Newton
				end
				jointLocation(j) = -1;
				LTE(j) = -1;
				fprintf('\t solving joint profile at station %g \n',station(stationRanges(j)))
                %update V2022-03-09 -> Get a rough estimate of the distance
                %between the tsd wheel and the joint
                %The distance between the current measurement and the joint's location should be a rough estimate of where the joint is at [+/- 10-15cm]
                estimateC = jointPosition(i) - station(stationRanges(j));  
                estimateC = max(estimateC-0.15,0):0.01:estimateC+0.15;
				
                %% solving the joint
                %Update V 2022-03-23 -> Use back-calc based on deflection velocity! May need vy and local vx record too!
                localvx = vx(stationRanges(j));
                
                %% IMPORTANT: VYcomes in mm/s, and vx in m/s. PASS BOTH IN M/sec
                 %Update v 2022-03-09 -> Pass an estimate of the joint location [variable estimateC]
                 verborragia = 1;
                [jointLocation(j),LTE(j),SSEfinal(j),~,~,~] =  backCalc_joint_BruteForce_0420(TSDpoints,localTSD_pulse(j,:)./1e3,localvx,estimateC,localThck,localK,localG,localE,nu,localLoad,pressure,verborragia);
				 
			end
			
			%% Update v 2022-03-04.
            %Do the export of the results now. I don't want to deal with the structure (that can be messy to export in few steps)
			%exportFilename stated at the beginning!
			exportSheet = sprintf('joint at station %g',jointPosition(i));
			%export k, E, G
			export = xlswrite(exportFilename,[localK;localE;localG],exportSheet,'e7:e9');
			%export h, nu
			export = xlswrite(exportFilename,[nu;localThck],exportSheet,'d13:d14');
			%export joint station, lat. and long.
			export = xlswrite(exportFilename,jointPosition(i),exportSheet,'f4');
			%         export = xlswrite(exportFilename,jointLatLong,exportSheet,'f2:f3');
			
			%export the c, LTE back-calc results upon approach
            %update v 2022-03-09 -> Export the final SSE
			exportVariable = [station(stationRanges),jointLocation,LTE,SSEfinal];
			export = xlswrite(exportFilename,exportVariable,exportSheet,'h6:k35');
            
            toc
		end  %end for i = 1 : nnn
		%% before closing the while, check if doing one more joint.		
		selectJoint = 0;
    end
save guardaTodo5cm_Loop124_v0420.mat
end 
disp('....all completed')
toc
%% RESTORE THE SINGULAR MATRIX WARNING
warning(warnStruct);