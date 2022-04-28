function [yd1,rT,WCT] = TSD_waveletDecomposition_WEIGHTED(X1,ydstart,f,LAMBDA,WeightWavelet,WeightDirac)
%function [yd1,rT,WCT] = TSD_waveletDecomposition(X1,ydStart,f,T)
%function to decompose a given 'ydStart' signal into a continuous (wavelet based) component yd, and a discontinuous component rT 
%and filter by thresholding with a threshold T.
%PROCEED BY BATCH GRADIENT DESCENT [SEE PAPER 1]
%input
%   ydStart = RAW SIGNAL [STRETCHED TO 2^P LENGTH]
%   f       = wavelet filter (provided from the TSD_denoisingJoints function)
%   LAMBDA  = threshold value [a.k.a. lambda]
%   WeightWavelet= diagonal matrix (length of ydStart) with the weights for
%   each component of the wavelet dictionary
%   WeightDirac  = diagonal matrix (length of ydStart) with the weights for
%   each component of the Dirac dictionary.
% 
%output
%   yd1     = continuous (wavelet-based) component of X1
%   rT      = discontinuous (spike) component of X1
%   WCT     = filtered wavelet coefficients for yd1
%
%Release candidate v2022-05-01

%Don't use Matlab's A^-1, cause it's computer intensive. Since the weight
%matrices are diagonal,extract the diag.elements, invert them, and then
%rebuild
aux =diag(WeightWavelet);
%aux = 1./aux;
% WeightWaveletInv = aux.*eye(length(aux));
WeightedWaveletLambda = LAMBDA.*aux;
% 
aux = diag(WeightDirac);
%aux = 1./aux;
% WeightDiracInv = aux.*eye(length(aux));
WeightedDiracLambda = LAMBDA.*aux;

%% code begins
%run 100 iterations
yd = ydstart;
for i = 1:100
    r  = X1-yd;   %residual of RAW yn minus DENOISED yd
    %start by wavelet-transform 'r': apply the wavelet decomposition using the filter 'f'. 
    %WC is the wavelet coefficients
    WC = FWT_PO(r,1,f);
    %syntax means: decompose the stretched-out X1, up to the max extent 
    %(the min extent is the power of 2 that makes the length of the vector X1).
    %UPDATE V2020-08-19: ADD THE WEIGHING FACTORS INTO THE SOFT-SHRINKING
    %FUNCTION
    WCT = softShrinking(WC,WeightedWaveletLambda);  %filter the noise with the softShrinking rule
    %e) invert the wavelet decomposition
    yd1 = IWT_PO(WCT,1,f);
    %recompute r:
    r = X1 - yd1; %theoretically this should be a mere spikes function (because I'm substracting a "clean" recovered sinusoidal signal)
    %remove the unwanted spikes on this r
    rT = softShrinking(r,WeightedDiracLambda);
    %and now close the loop by redefining yd as the thresholded residual. ?Weird, ain't it?
    yd = rT;    
end

end %endfunction