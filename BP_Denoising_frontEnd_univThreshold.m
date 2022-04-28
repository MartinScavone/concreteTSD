%TSD De-noising procedure front-end code
%V0.1 2020-03-05
%This front-end script will aid to load the raw TSD data and make all
%necessary calls to functions that do the denoising procedure.
%
%Release candidate v2022-05-01

%%
tic
clc
clear variables
close all
addpath('./WAVELAB850')
addpath('./WAVELAB850/Orthogonal');  
addpath('./WAVELAB850/Utilities');  
% addpath('./export_fig');

disp('Loading TSD Data')

inputFilename = '0033_FHWA_EFL_PFS_1m_Deflection_only_6_13_19.xlsb';
inputSheet    = '0033_FHWA_EFL_PFS_1m_TSD';


latLongData = xlsread(inputFilename,inputSheet,'bj2:bo79257');  %%lat long data in decimal degrees
stationData = xlsread(inputFilename,inputSheet,'i2:j79257');    %% station in km
[~ ,roadID,~] = xlsread(inputFilename,inputSheet,'a2:a79257');   %%roadID is an alphanumeric field. STORED AS CELL ARRAY! this one combines the road id + the segm. id. Try as unique identifier...
[~,roadName,~] = xlsread(inputFilename,inputSheet,'e2:e79257'); %%roadName is an alphanumeric field. STORED AS CELL ARRAY!
goodMeasures= xlsread(inputFilename,inputSheet,'k2:y79257');    %% metric units.
%columns are:: SCI_200 (microns)	SCI_300 (microns)	SCI_SUBGRADE (microns)	D0 (microns)	D203 (microns)	D305 (microns)	D457 (microns)	D610 (microns)	D914 (microns)	D1219 (microns)	D1524 (microns)	D1829 (microns)

latFrom = latLongData(:,1);
latTo   = latLongData(:,3);
longFrom= latLongData(:,2);
longTo  = latLongData(:,4);
altFrom = latLongData(:,5);
altTo   = latLongData(:,6);

stationFrom = 1000*stationData(:,1);  %%get station information (from and to) for each measurement. Convert from km to m.
stationTo   = 1000*stationData(:,2);

metricD0 = goodMeasures(:,4);
metricD0 = -1*metricD0;

%% Get all unique roads in the TSD dataset - UPDATE V2019-10-30: MUST DISTINGUISH BY roadName+blockID [roadID and roadNAME by themselves would collect many segments together (and overlap signals)]
disp('National Mall Measurements: Analyzing all streets and segments separately')
roadIDList = unique(roadID);    %use this vector to locate the TSD points for each road.
numRoads = length(roadIDList);
fprintf('\t %g streets (+ blocks) were recognized \n',numRoads)

%store all the weak spot locations for all roadIDs in this struct. WIll unpack it later to map them all...
weakSpotsMap = struct('roadID',[],'roadName',[],'weakSpotLocations',[],'weakSpotSpac',[],'weakSpotLatFrom',[],'weakSpotLatTo',[],'weakSpotLongFrom',[],'weakSpotLongTo',[],'RawDo',[]); 


%% 3) WAVELET DENOISING OF THE TSD DATA - DECOMPOSITION INTO WAVELETS AND APPLICATION OF A SOFT-SHRINKING THRESHOLD
disp('TSD: Applying wavelet denoising with soft-shrinking threshold')
    
for k = 1:numRoads

   whichRoad = roadIDList(k);
   
   tsdIndices = find(strcmp(roadID,whichRoad));  %this tweak allows to use find to compare cell arrays. Source: https://www.mathworks.com/matlabcentral/answers/84242-find-in-a-cell-array
   %get the name of the road/lane under analysis
   thatRoadName = roadName(tsdIndices);
   thatRoadName = thatRoadName(1);            
   
   %get the TSD data for that road
   fprintf(' \t now processing: %s, segment %s \n',string(thatRoadName),string(whichRoad));
   refTSD = metricD0(tsdIndices);
   station = stationFrom(tsdIndices);
   latStart = latFrom(tsdIndices);
   latEnd   = latTo(tsdIndices);
   longStart= longFrom(tsdIndices);
   longEnd  = longTo(tsdIndices);

 
    %%sanity control: the refTSD vectors may contain NaNs.
    %DON'T DO::: force replace with zeros...  <- THIS DROPS STDEV TO ZERO %(WHEN THE ACTUAL DATA HAS MORE VARIABILITY!)
    %DO::::::::: locate the nans, remove them altogether and forget about them!  (Samer's advice 2019-11-04)
    
    %use this workaround: https://www.mathworks.com/matlabcentral/answers/164316-select-everything-not-returned-by-index    
    theseAreNans = find(isnan(refTSD)); %this returns a vector with the positions of NaN's inside refTSD
    NansIndex = isnan(refTSD);  %this logical array (different from the one above) helps then get the non-NaN values easily. [it's a 0/1 vector as long as refTSD in which the 1's indicate the NaNs]

    
    %sanity control 02: it may happen that all the refTSD points are or blanks [there's actually a case out there].
    %Kill the iteration if such is the case (fill the weakspotMap(k) and do nothing else....
    if length(theseAreNans) == length(refTSD)
        weakSpotsMap(k).roadID = whichRoad;
        weakSpotsMap(k).roadName = thatRoadName;
        weakSpotsMap(k).weakSpotSpac      = 0;
        weakSpotsMap(k).weakSpotLocations = 0;
        weakSpotsMap(k).weakSpotLatFrom   = 0;
        weakSpotsMap(k).weakSpotLongFrom  = 0;
        weakSpotsMap(k).weakSpotLatTo     = 0;
        weakSpotsMap(k).weakSpotLongTo    = 0;
        weakSpotsMap(k).rawDo             = 0;
        continue
    end
    
    %case this not occurs, remove the theseAreNans and keep on with the    %calculations.
    refTSD = refTSD(~NansIndex);
    station = station(~NansIndex);
    latStart = latStart(~NansIndex);
    latEnd   = latEnd(~NansIndex);
    longStart= longStart(~NansIndex);
    longEnd  = longEnd(~NansIndex);
    
    n = length(refTSD);   %length of the sanitized series.
    
    %% Iterative Basis pursuit decomposition + soft thresholding procedure starts.... 
%   pass to the TSD_denoisingJoints function!
    callType = 'default';   %No other callType allowable!

 %%Pass data for the plot: T vs the steinError and locate the optimum T value (that who %minimizes the SURE)
    %the plot is made at the TSD_denoisingJoints function. 
    SUREPlotNumber = -1;
    SUREPlotName   = 'pototo';
%     SUREPlotName   = sprintf('Denoised TSD - %s, section %s : optimum t value location',string(thatRoadName), string(whichRoad));
    
    [denoisedTSD,~,ySpikes,optLambda,optSURE,SUREPlotHandle] = TSD_denoisingJoints_UnivThrshld(refTSD,callType,SUREPlotNumber,SUREPlotName);

    %% 5) Figure 100 series with the final denoised signal obtained with the optimum Threshold value.

    figure(100+k)
    plotName = sprintf('Denoised TSD - %s, section %s : Noisy signal vs. denoised signal',string(thatRoadName), string(whichRoad));
    subplot(2,1,1)
    plot (station, refTSD, 'color',[0.73 0.83 0.96])
    set(gca, 'FontName', 'Arial')
    set(gca,'defaultAxesFontSize',16)
    grid on
    titleString = sprintf('TSD signal for road %s, section %s',string(thatRoadName), string(whichRoad));
    title(titleString)
    xlabel ('station [m]')
    ylabel ('D_0 [\mu m]')
    hold on
    plot(station,denoisedTSD,'color',[0.93 0.69 0.13]) %orange color
    legend('noisy','Denoised BP','Denoised Rew L1')
    hold off
    
    subplot(2,1,2)
    plot(station,ySpikes,'bo','linewidth',1,'markersize',6)
    set(gca, 'FontName', 'Arial')
    set(gca,'defaultAxesFontSize',16)
    grid on
    title('Recovered Dirac component. Universal threshold')
    legend('BP recovered')
    
    %send the plot to a pdf file with export_fig.  
%     export_fig NatlMallTSD_univThreshold -pdf -append   
    
    %% 6) Add this final analysis to detect how many spikes were detected
    %and how many of those are real.

    id2 = find(ySpikes ~=0);
    fprintf('\t \t %g spikes were detected with the denoising \n',length(id2));

    %get the spacing between weak spots and see if they correspond with the
    %actual joint spacing [12.5m, Katicha et al., 2013]
    weakSpotsPositions = find(ySpikes>0);  %   %%note that vector position does match the station and distance betwee nweekspots because each tsd measurement is spaced 1.00m
    weakSpotSpacing = weakSpotsPositions(2:end)-weakSpotsPositions(1:end-1);
    meanJointSpacing = mode(weakSpotSpacing);  %note: prefer mode and/or median to the arithmetic mean of weakSpotSpacing because the mean may be stretched to a too large value by long segments w/o spots
    
    %save results for mapping...
    weakSpotsMap(k).roadID = whichRoad;
    weakSpotsMap(k).roadName = thatRoadName;
    weakSpotsMap(k).weakSpotSpac      = meanJointSpacing;
    weakSpotsMap(k).weakSpotLocations = station(weakSpotsPositions);
    weakSpotsMap(k).weakSpotLatFrom   = latStart(weakSpotsPositions);
    weakSpotsMap(k).weakSpotLongFrom  = longStart(weakSpotsPositions);
    weakSpotsMap(k).weakSpotLatTo     = latEnd(weakSpotsPositions);
    weakSpotsMap(k).weakSpotLongTo    = longEnd(weakSpotsPositions);
    weakSpotsMap(k).RawDo             = refTSD(weakSpotsPositions);

    %%
%     figure(400+k)
%     plotName = sprintf('Denoised TSD - %s, section %s : Spacing between weak points',string(thatRoadName), string(whichRoad));
%     set(gcf,'Name',plotName);
%     hist(weakSpotSpacing,20);
%     title('Distribution of weak spot separations')
%     xlabel('weak spot separation [m]')
%     ylabel('count')
     
     %send the plot to a pdf file with export_fig.
     %export_fig NatlMallTSD -pdf -append
     fprintf('\t Road %g of %g completed \n',k,numRoads)

end  %%end-for-length(roadnames)

%% cleanup...

clear weakSpotsPositions
clear weakSpotSpacing
clear minSteinPos
clear steinError
clear denoisedTSD
clear peaksMaps

clear refTSD station latEnd latStart longEnd longStart tsdIndices

%% export the weakSpotsMap to Excel for easy mapping
%use the already started file
exportFilename = 'NatlMall_weakSpotChart_univV2022-0501.xlsx';
disp('exporting Weak Spots Map to Excel')
startRow = 3;
%end row = startRow + length-1
%and new start row = endrow + 1
for i = 1:length(weakSpotsMap)
   fprintf(' \t Progress %g of %g \n',i,length(weakSpotsMap))
   %export all struct contents one by one
   len = length(weakSpotsMap(i).weakSpotLocations);
   endRow = startRow+len-1;
   exportVariable=[weakSpotsMap(i).weakSpotLocations weakSpotsMap(i).weakSpotSpac*ones(len,1) weakSpotsMap(i).weakSpotLatFrom weakSpotsMap(i).weakSpotLatTo weakSpotsMap(i).weakSpotLongFrom weakSpotsMap(i).weakSpotLongTo weakSpotsMap(i).RawDo];
   exportRange = strcat('d',string(startRow),':j',string(endRow));
   %export all numerics
   uu = xlswrite(exportFilename,exportVariable,1,exportRange);
   %export cell contents.
   exportRange = strcat('b',string(startRow),':c',string(endRow));
   uu = xlswrite(exportFilename,[weakSpotsMap(i).roadID weakSpotsMap(i).roadName],1,exportRange);
   
    %done exporting. set the startRow forthe next export.
    startRow = endRow+1;
end
disp('....all completed')
toc