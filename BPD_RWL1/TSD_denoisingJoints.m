function [yCleanExp,WCTsExp,peaksMapsExp,optLambda,steinErrorExp,SUREPlotHandle] = TSD_denoisingJoints(yRaw,callType,doThePlot,SUREPlotNumber,SUREPlotName)
%function [yCleanExp,WCTsExp,peaksMapsExp,optLambda,steinErrorExp,SUREPlotHandle] = TSD_denoisingJoints(yRaw,callType,doThePlot,SUREPlotNumber,SUREPlotName)
%
%Remove noise from the TSD measurements by basis pursuit and soft shrinking
%thresholding. Main function. 
%This function receives a noisy TSD measurement and a type of wavelet
%filter to select [CONSULT WAVELAB HELP] and applies the denoising
%procedure [Katicha et al., 2013] over a different range of possible
%threshold values (lambdaRange)
%
%Input:
%   yRaw:           Raw TSD measurement [mandatory!]
%   callType        Either 'default' [use default wavelet (Symmlet 8th degr.) and threshold option], or 'custom'[ enter your own]
%   SUREPlotNumber  Number to identify the figure with the SURE(lambda) plot
%   SUREPlotName    Name to give to that figure
%   doThePlot       Boolean to tell the code if to do the SURE plot (1) or
%                   not (0)
%
%Output:
%
%Dependencies:
%This function depends on the following:
%   WAVELAB library - addpath'd from the FrontEnd
%   TSD_waveletDecomposition function%
%
%Candidate release version v2022-05-01
%
	
%% Preprocessing: set defaults.

%prepare the denoising procedure
if strcmp(callType,'default')
    %use defaults
    %a) compute signal standard deviation [refer to Katicha et al., 2015], %equation 11
    %Update v2022-02-12 - MAD formula for stdev(diff(y)) and built-in MAD function
    auxRaw = diff(yRaw);
    sigma = mad(auxRaw,1); 
    sigma = 1.4826.*sigma;
    sigma = sigma./sqrt(2);
   
    %B) default wavelet type 
    waveletType = 'Symmlet';
    waveletOrder = 8;
 
 
else
    error('TSD_denoisingJoints error: invalid call type')
end

%% - Preprocessing: need to expand yRaw to 2^n length
%normalize yRaw [use either user-provided sigma or default sigma estimation]
yRaw = yRaw/sigma;
n = length(yRaw);

%compute lambda range for the denoising (must be done over the unfolded
%vector. It is user-provided if callType = 'custom'
if strcmp(callType,'default')
    lambdaRange = 0.0:0.05:sqrt(2*log(length(yRaw)));
    lambdaRange(1) = 0.01;
else
    %do Nothing, lambdaRange is passed as input.
end

%unfold yRaw to 2^n size
[targetSize,~] = powerOf2(n);
X1 = unfoldVector(yRaw, targetSize);   

%define the wavelet filter
f = MakeONFilter(waveletType,waveletOrder);
%default: f = MakeONFilter('Symmlet',8);


%% - Processing: initialize output variables

%store all values of the denoised TSD dataset with different threshold and
%the estimated error (using Stein's SURE estimate [Samer's suggestion])

X1T = zeros(length(X1),length(lambdaRange));            %< unfolded filtered tsd data
yd1 = zeros(length(X1),length(lambdaRange));            %< unfolded wavelet component
rT  = zeros(length(X1),length(lambdaRange));            %< unfolded discontinuous peak component
WCT = zeros(length(X1),length(lambdaRange));            %< sanitized wavelet coefficients for the continuous component.

yClean = zeros(n,length(lambdaRange));                  %< folded-back filtered tsd data
steinError = zeros(size(lambdaRange));                  %< stein's SURE error estimate for each value of T  
peaksMaps  = zeros(n,length(lambdaRange));              %< folded-back spikes' component

%flip over the lambda vector (so that it starts with the largest threshold
lambdaRange = flip(lambdaRange);

%% - Processing: Iterative Gradient Descent and SURE evaluation for each value of the threshold Lambda

for tt = 1:length(lambdaRange)
    T = lambdaRange(tt);
    if tt == 1 
    %initialize the gradient descent. Start with a vector of zeros in the
    %first iteration over T, or use the clear signal for the T(i+1) when
    %doing T(i) 
        ydStart = zeros(size(X1));
    else
        ydStart = X1T(:,tt-1);
    end
    %Decompose as wavelets and do the Batch Gradient Descent.
    [yd1(:,tt),rT(:,tt),WCT(:,tt)] = TSD_waveletDecomposition(X1,ydStart,f,T);  
    
    %The denoised signal is the denoised spikes (rT) plus the denoised sinusoidal component
    X1T(:,tt) = rT(:,tt) + yd1(:,tt);     

    %Compute the SURE (Stein's unbiased risk estimate) - stored as steinError.
    steinError(tt) = TSD_SURE(X1, X1T(:,tt),WCT(:,tt),rT(:,tt));

    %Roll back the streched-out denoised signal X1T and muitiply by sigma to restore it to TSD "units"
    yClean(:,tt) = foldBackVector(X1T(:,tt),n);
    yClean(:,tt) = yClean(:,tt)*sigma;
    peaksMaps(:,tt)   = foldBackVector(rT(:,tt),n);
    peaksMaps(:,tt)   = peaksMaps(:,tt)*sigma;   %Retrieve the discontinuous component, fold back, and re-set to scale.
end


%% Locate the optimum lambda value - position of the min(steinError).
minSteinPos = find(steinError == min(steinError));
if ~isempty(minSteinPos)
    minSteinPos = minSteinPos(1);
else
    minSteinPos = 1;  %force some value in case a flat signal appears (there is one such case)
end
   
optLambda = lambdaRange(minSteinPos);

%%Export: 
steinErrorExp = steinError(minSteinPos);
yCleanExp = yClean(:,minSteinPos);
peaksMapsExp = peaksMaps(:,minSteinPos);  
WCTsExp = WCT(:,minSteinPos);

%% Plot result
%create the figure of SURE versus lambda and pass the handle to the
%frontEnd
if doThePlot
    SUREPlotHandle = figure(SUREPlotNumber);
    set(SUREPlotHandle,'name',SUREPlotName)
    plot(lambdaRange,steinError,'b')
    grid on
    title('TSD noise reduction: SSE prediction with SURE')
    xlabel('Soft shrinking threshold' )
    ylabel('SURE estimate')

    hold on
    plot(optLambda,steinErrorExp,'r+','markersize',12)
    legend('SURE estimate','Optimum SURE')
    hold off
else
    SUREPlotHandle = -1;
end
    
end  %end-function
