function SURE = TSD_SURE(X1, X1T,WCT,rT)
%function SURE = TSD_SURE(X1, X1T,WCT,rT)
%
%Compute the SURE estimate of the BPD procedure for a given series X1, and
%its filtered components X1T, and the vectors of coefficeints WCT, rT
%
%Release candidate v 2022-05-01 

%First, Count how many WCT(:,k) coefficients are non-zero
%%THE STEIN ERROR TERM MUST INCLUDE THE COUNT OF NON-ZERO SHRUNK COEFFICIENTS FOR THE TWO SOFT-SHRINKING STEPS (SINUSOID + SPIKES)!
nonZeroWCTcount1 = length(find(WCT~=0));  
nonZeroWCTcount2 = length(find(rT ~=0));
%
%%IMPORTANT: NEED TO COMPUTE TEH STEIN ERROR ESTIMATE WITH THE X1 AND X1T SERIES (THE FOLDBACK TAKES SOME OF THE OPTIMIZATION AWAY AND MAY RETURN A FAULTY DENOISED SIGNAL!     
SURE = -length(X1T)+sum((X1T-X1).^2)+2*nonZeroWCTcount1+2*nonZeroWCTcount2;   


end