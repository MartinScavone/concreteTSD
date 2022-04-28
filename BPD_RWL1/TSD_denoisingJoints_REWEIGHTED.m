function [yClean,waveletMaps,spikeMaps,steinError,optiLambda] = TSD_denoisingJoints_REWEIGHTED(yRaw,lambdaRange,iterations,weightsEpsilon,optResultsOnly)
%function [yClean,waveletMaps,spikeMaps,steinError,optimLambda] = TSD_denoisingJoints_REWEIGHTED(yRaw,lambdaRange,iterations,weightsEpsilon,optResultsOnly)
%
%Remove noise from the TSD measurements by Reweighted L1 minimization.
%This function receives a noisy TSD measurement and proceeds iteratively
%as prompted in the paper on Reweighted L1 Min. by Candés et al 2008
%
%Input:
%   yRaw:           Raw TSD measurement [mandatory!]
%   lambdaRange:    Give a vector of values for the penalty parameter Lambda. Or pass -1 to seek the optimum-fitting lambda (call's the TSD_denoisingJoints function)
%   iterations:     Number of loops of the reweighted L1 minimization
%   weightsEpsilon: Stability parameter (eps) -> See candès et al., 2008
%   optResultsOnly: Boolean value to tell if the code shall return
%   yClean/and components for all values of lambdaRange (=0) or only for
%   the SURE-optimizing case (=1)
%Output:
%   yClean:         Denoised signal
%   waveletMaps:    Coefficients for the continuous component onto the
%   wavelet signal space
%   spikeMaps:      Coefficients for the discontinuous component onto the
%   Dirac signal space
%   steinError:     Values of SURE [SSE estimation] for each value of
%   lambdaRange
%   IMPORTANT: If optResultsOnly = 1, the output is cut to the end results
%   for the SURE optimizng case only
%
%Dependencies:
%This function depends on the following:
%   WAVELAB library - addpath'd from the FrontEnd
%   TSD_waveletDecomposition function
%
%Release candidate V2022-05-01
%
%V1.2 2022-02-18:
%Stability update - when exporting, if more than one value of lambda
%minimizes SURE, report one only. It's supposed never to happen, but it
%somehow happens on some special cases.
%
%V1.1 2021-05-10:
%-added the last stage of correction using the improved SURE formula to
%select the optimal lambda value
%(Mazumder et al., 2011), it appears in Paper 2.


%% Preprocessing: set defaults.
%Update v2022-02-12 - MAD formula corrected by Samer (He was actually right)
%-- use built-in MAD function
auxRaw = diff(yRaw);
sigma = mad(auxRaw,1); 
sigma = 1.4826.*sigma;
sigma = sigma./sqrt(2);

waveletType = 'Symmlet';
waveletOrder = 8;

% weightsEpsilon = 2;      
%% - Preprocessing: need to expand yRaw to 2^n length
%normalize yRaw [use either user-provided sigma or default sigma estimation]
yRaw = yRaw/sigma;
n = length(yRaw);

%unfold yRaw to 2^n size
[targetSize,~] = powerOf2(n);
X1 = unfoldVector(yRaw, targetSize);   

%% define the wavelet filter
f = MakeONFilter(waveletType,waveletOrder);
%default: f = MakeONFilter('Symmlet',8);

%% - Processing: initialize output variables
%store all values of the denoised TSD dataset with different threshold and
%the estimated error (using Stein's SURE estimate [Samer's suggestion])

X1T = zeros(length(X1),iterations,length(lambdaRange));            %< unfolded filtered tsd data
yd1 = zeros(length(X1),iterations,length(lambdaRange));            %< unfolded wavelet component
rT  = zeros(length(X1),iterations,length(lambdaRange));            %< unfolded discontinuous peak component
WCT = zeros(length(X1),iterations,length(lambdaRange));            %< sanitized wavelet coefficients for the continuous component.

steinError = zeros(1,iterations,length(lambdaRange));
optiLambda = lambdaRange;

yClean = zeros(n,iterations,length(lambdaRange));                  %< folded-back filtered tsd data
spikeMaps = zeros(n,iterations,length(lambdaRange));               %< folded-back peaks
waveletMaps=zeros(n,iterations,length(lambdaRange));               %< folded back wavelet component?   

%initialize the weight vectors
weightsWavelets = zeros(length(X1),iterations,length(lambdaRange));
weightsDiracs   = zeros(length(X1),iterations,length(lambdaRange));
%fill up the first iteration for all lambdaRanges
weightsWavelets(:,1,:) = ones(length(X1),1,length(lambdaRange));
weightsDiracs(:,1,:)   = ones(length(X1),1,length(lambdaRange));

%% NEW: INITIALIZE WEIGHT MATRICES - Prepare separate matrices for 
    %% - Processing: Iterative Gradient Descent and SURE evaluation for each value of the threshold Lambda    
for tt = 1:length(lambdaRange)
    lambda = lambdaRange(tt);
    if tt == 1 
    %initialize the gradient descent. Start with a vector of zeros in the
    %first iteration over T, or use the clear signal for the T(i+1) when
    %doing T(i) 
        ydStart = zeros(size(X1));
    else
        ydStart = X1T(:,iterations,tt-1);
    end
    
    for k = 1:iterations
        %the first iteration is " plain BPD" - initialize the diagonal
        %matrix of weightsas ID
        WWCT = eye(length(X1));
        Wspk = eye(length(X1));
        %fill up the weights onto the diagonals of the W matrices
        WWCT = reshape(weightsWavelets(:,k,tt),[length(X1),1]).*WWCT;
        Wspk = reshape(weightsDiracs(:,k,tt),[length(X1),1])  .*Wspk;        

        %Decompose as wavelets and do the Batch Gradient Descent.
        [yd1(:,k,tt),rT(:,k,tt),WCT(:,k,tt)] = TSD_waveletDecomposition_WEIGHTED(X1,ydStart,f,lambda,WWCT,Wspk);  

        %The denoised signal is the denoised spikes (rT) plus the denoised wavelet component
        X1T(:,k,tt) = rT(:,k,tt) + yd1(:,k,tt);     
        
        %Compute the SURE (Stein's unbiased risk estimate) - stored as steinError.
        steinError(1,k,tt) = TSD_SURE(X1, X1T(:,k,tt),WCT(:,k,tt),rT(:,k,tt));

        %Roll back the streched-out denoised signal X1T and muitiply by sigma to restore it to TSD "units"
        yClean(:,k,tt) = foldBackVector(X1T(:,k,tt),n);
        yClean(:,k,tt) = yClean(:,k,tt).*sigma;
        spikeMaps(:,k,tt)   = foldBackVector(rT(:,k,tt),n);
        spikeMaps(:,k,tt)   = spikeMaps(:,k,tt).*sigma;   %Retrieve the discontinuous component, fold back, and re-set to scale.
        waveletMaps(:,k,tt) = foldBackVector(yd1(:,k,tt),n);
        waveletMaps(:,k,tt) = waveletMaps(:,k,tt).*sigma;
        
        %% finally, prepare the coefficients for the next iteration
        %avoid this step if k == iterations
        if k< iterations
           weightsWavelets(:,k+1,tt) = 1./(abs(WCT(:,k,tt)) + weightsEpsilon);
           weightsDiracs(:,k+1,tt)   = 1./(abs(rT(:,k,tt))  + weightsEpsilon); 
        end
    end %end loop for k (given value of lambda)
    
end %end loop for all values of lambda

%% update V. 2021-05-10
% Add here the improved SURE value for all lambdas 
% This is to run after all iterations completed

steinUpdated = zeros(1,length(lambdaRange));

for tt = 1:length(lambdaRange)
    %compute the updated SURE error estimate using the fit from after the
    %RWL1 (X1T(iterations)) and the count of non-zero components from the BP fit [first
    %iteration] -> See eqn. "11?" in paper 2.
    steinUpdated(tt) = TSD_SURE(X1, X1T(:,iterations,tt),WCT(:,1,tt),rT(:,1,tt));
end

%replace the steinError output with steinUpdated one...
steinError = steinUpdated;

%also, make the output signals correspond to the steinUpdated-minimizing
%case, and the last iteration of RWL1 (I don't mind about the results
%half-way-through)
if optResultsOnly
    optLambdaPos = steinUpdated==min(steinUpdated);
    if sum(optLambdaPos)>1  %use sum to check if more than one 1, cause optLambdaPos is a BOolean vector the size of steinUpdated.
        aux = find(optLambdaPos);
        aux = aux(1);
        bux = zeros(size(optLambdaPos));
        bux(aux) = 1;        
        optLambdaPos = boolean(bux);
    end    
    optiLambda  = lambdaRange(optLambdaPos);
    yClean      = squeeze(yClean(:,end,optLambdaPos));
    spikeMaps   = squeeze(spikeMaps(:,end,optLambdaPos));
    waveletMaps = squeeze(waveletMaps(:,end,optLambdaPos));   
    steinError = min(steinError);
end

end    %endfunction
