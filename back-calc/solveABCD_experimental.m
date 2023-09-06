function [A,B,C,D] = solveABCD_experimental(yDomain,sDomain,p,a,b,k,g,l,nu,delta,c,verboseness)
%function [A,B,C,D] = solveABCD_experimental(yDomain,sDomain,p,a,b,k,g,l,nu,delta,c,verboseness)
%
%Auxiliary function that solves the boundary conditions for the slab with a joint problem 
%IMPLEMENTATION WITH SLAB THEORY - VAN CAUW. CHAPTER 15
%V2.0 experimental 2022-03-10
%   Re-computed the derivatives in a tidier way (it's causing countless
%   approximation errors and NaNs when computing)

%% PREPARATION OF THE SYSTEM OF EQUATIONS - SHOULD BE 4 X 4 x length of the S domain
%State it as M X = N
%Prepare the matrices M and N by row - each row represents one equation
sss = length(sDomain); 
M = zeros(4,4,sss);
N = zeros(4,1,sss);

A = zeros(sss,length(yDomain));
B = zeros(sss,length(yDomain));
C = zeros(sss,length(yDomain));
D = zeros(sss,length(yDomain));

cPlus =1.01.*c;
cMinus=0.99.*c;
c2plus=1.02.*c;
c2minus=0.98.*c;

deltaX = cPlus-c;
if length(yDomain)>1
    deltaY = 0.3.*(yDomain(2)-yDomain(1));
    yDomainPlus = yDomain+0.3.*deltaY;
    yDomainMinus= yDomain-0.3.*deltaY;
else
    deltaY = 0.01;
    yDomainPlus = yDomain + deltaY;
    yDomainMinus = yDomain - deltaY;
end

for j = 1:length(yDomain)
    if verboseness
        fprintf(' \t \t SolveABCD: Computing boundary conditions %g of %g \n',j,length(yDomain))
    end        
    for i = 1:sss
        s = sDomain(i);
        %I may need these quantities for equation 1
        wX = compute_w_term(c,yDomain(j),s,p,a,b,k,g,l);
        [waX,wbX] = compute_wawb_term(c,yDomain(j),s,p,b,k,g,l);
        [wcX,wdX] = compute_wcwd_term(c,yDomain(j),s,p,b,k,g,l);
        
        %% Prepare equation 1 - LTE (eq. 15.35)
        M(1,:,i) = [delta*waX delta*wbX -wcX -wdX];
        N(1,:,i) = -delta*wX;

       %I may need these guys for equations 2-4 (note: must calculate over all
       %yDomain so that I can derive afterwards)
       wXplus = compute_w_term(cPlus,yDomain(j),s,p,a,b,k,g,l);
       [waXplus,wbXplus] = compute_wawb_term(cPlus,yDomain(j),s,p,b,k,g,l);
       [wcXplus,wdXplus] = compute_wcwd_term(cPlus,yDomain(j),s,p,b,k,g,l);
       
       wXminus = compute_w_term(cMinus,yDomain(j),s,p,a,b,k,g,l);
       [waXminus,wbXminus] = compute_wawb_term(cMinus,yDomain(j),s,p,b,k,g,l);
       [wcXminus,wdXminus] = compute_wcwd_term(cMinus,yDomain(j),s,p,b,k,g,l);
       
       wYplus = compute_w_term(c,yDomainPlus(j),s,p,a,b,k,g,l);
       [waYplus,wbYplus] = compute_wawb_term(c,yDomainPlus(j),s,p,b,k,g,l);
       [wcYplus,wdYplus] = compute_wcwd_term(c,yDomainPlus(j),s,p,b,k,g,l);
       
       wYminus = compute_w_term(c,yDomainMinus(j),s,p,a,b,k,g,l);
       [waYminus,wbYminus] = compute_wawb_term(c,yDomainMinus(j),s,p,b,k,g,l);
       [wcYminus,wdYminus] = compute_wcwd_term(c,yDomainMinus(j),s,p,b,k,g,l);
       
       %these guys Id need for the 3rd derivative over x
       wX2plus = compute_w_term(c2plus,yDomain(j),s,p,a,b,k,g,l);
       [waX2plus,wbX2plus] = compute_wawb_term(c2plus,yDomain(j),s,p,b,k,g,l);
       [wcX2plus,wdX2plus] = compute_wcwd_term(c2plus,yDomain(j),s,p,b,k,g,l);
       
       wX2minus = compute_w_term(c2minus,yDomain(j),s,p,a,b,k,g,l);
       [waX2minus,wbX2minus] = compute_wawb_term(c2minus,yDomain(j),s,p,b,k,g,l);
       [wcX2minus,wdX2minus] = compute_wcwd_term(c2minus,yDomain(j),s,p,b,k,g,l);
       
       %these guys I need for the d3w/dxdy2
       wXplusYplus = compute_w_term(cPlus,yDomainPlus(j),s,p,a,b,k,g,l);
       [waXplusYplus,wbXplusYplus] = compute_wawb_term(cPlus,yDomainPlus(j),s,p,b,k,g,l);
       [wcXplusYplus,wdXplusYplus] = compute_wcwd_term(cPlus,yDomainPlus(j),s,p,b,k,g,l);
       
       wXplusYminus = compute_w_term(cPlus,yDomainMinus(j),s,p,a,b,k,g,l);
       [waXplusYminus,wbXplusYminus] = compute_wawb_term(cPlus,yDomainMinus(j),s,p,b,k,g,l);
       [wcXplusYminus,wdXplusYminus] = compute_wcwd_term(cPlus,yDomainMinus(j),s,p,b,k,g,l);
       
       wXminusYplus = compute_w_term(cMinus,yDomainPlus(j),s,p,a,b,k,g,l);
       [waXminusYplus,wbXminusYplus] = compute_wawb_term(cMinus,yDomainPlus(j),s,p,b,k,g,l);
       [wcXminusYplus,wdXminusYplus] = compute_wcwd_term(cMinus,yDomainPlus(j),s,p,b,k,g,l);
       
       wXminusYminus = compute_w_term(cMinus,yDomainMinus(j),s,p,a,b,k,g,l);
       [waXminusYminus,wbXminusYminus] = compute_wawb_term(cMinus,yDomainMinus(j),s,p,b,k,g,l);
       [wcXminusYminus,wdXminusYminus] = compute_wcwd_term(cMinus,yDomainMinus(j),s,p,b,k,g,l);
       
      %% derive w(x) three times
       
      %update v 2022-03-12 : re-write the derivative as increment over plus - minus / 2deltax
      %approximate it directly to x == c
      dwdx = (wXplus - wXminus)./(2.*deltaX);  %%this is the derivative AT x == c AND ydomain(j)!
%       dwdy = (wYplus - wYminus)./(2.*deltaY);  %%this is the derivative AT x == c AND ydomain(j)!
      
      %update v2022-03-12 -> do a second-order-central derivative for d2w/dx2 and d2w/dy2. 
      %Directly estimated at x == c and  all yDomain
      d2wdx2 = (wXplus + wXminus - 2.*wX)./(deltaX.^2);  %these results are valid for x == c and AND ydomain(j)!
      d2wdy2 = (wYplus + wYminus - 2.*wX)./(deltaY.^2);
      
      %update v2022-03-12: solve the third derivatives as single derivatives of the d2wdy2
      d2wdx2plus = (wX2plus - 2.*wXplus + wX)./(deltaX.^2);
      d2wdx2minus= (wX - 2.*wXminus + wX2minus)./(deltaX.^2);
      d2wdy2plus = (wXplusYplus - 2.*wXplus + wXplusYminus)./(deltaY.^2);
      d2wdy2minus= (wXminusYplus - 2.*wXminus + wXminusYminus)./(deltaY.^2);
            
      d3wdx3 = (d2wdx2plus - d2wdx2minus)./(2.*deltaX);      
      d3wdxdy2=(d2wdy2plus - d2wdy2minus)./(2.*deltaX);
      
      %----------------
      %same story for deriving wa(x) 3 times -> UPDATED V2022-03-12
      dwadx = (waXplus - waXminus)./(2.*deltaX);
%       dwady = (waYplus - waYminus)./(2.*deltaY);	
      
      d2wadx2 = (waXplus + waXminus - 2.*waX)./(deltaX.^2);
      d2wady2 = (waYplus + waYminus - 2.*waX)./(deltaY.^2);
      
      d2wAdx2plus = (waX2plus - 2.*waXplus + waX)./(deltaX.^2);
      d2wAdx2minus= (waX - 2.*waXminus + waX2minus)./(deltaX.^2);      
      d2wAdy2plus = (waXplusYplus - 2.*waXplus + waXplusYminus)./(deltaY.^2);
      d2wAdy2minus= (waXminusYplus - 2.*waXminus + waXminusYminus)./(deltaY.^2);
      
      d3wadx3 = (d2wAdx2plus - d2wAdx2minus)./(2.*deltaX);      
      d3wadxdy2=(d2wAdy2plus - d2wAdy2minus)./(2.*deltaX);
      
      %----------------
      %same story for deriving wb(x) 3 times -> UPDATED V2022-03-12
      dwbdx = (wbXplus - wbXminus)./(2.*deltaX);
%       dwbdy = (wbYplus - wbYminus)./(2.*deltaY);	
      
      d2wbdx2 = (wbXplus + wbXminus - 2.*wbX)./(deltaX.^2);
      d2wbdy2 = (wbYplus + wbYminus - 2.*wbX)./(deltaY.^2);
      
      d2wBdx2plus = (wbX2plus - 2.*wbXplus + wbX)./(deltaX.^2);
      d2wBdx2minus= (wbX - 2.*wbXminus + wbX2minus)./(deltaX.^2);      
      d2wBdy2plus = (wbXplusYplus - 2.*wbXplus + wbXplusYminus)./(deltaY.^2);
      d2wBdy2minus= (wbXminusYplus - 2.*wbXminus + wbXminusYminus)./(deltaY.^2);
      
      d3wbdx3 = (d2wBdx2plus - d2wBdx2minus)./(2.*deltaX);      
      d3wbdxdy2=(d2wBdy2plus - d2wBdy2minus)./(2.*deltaX);
      
      %----------------
      %same story for deriving wc(x) 3 times -> UPDATED V2022-03-12
      dwcdx = (wcXplus - wcXminus)./(2.*deltaX);
%       dwcdy = (wcYplus - wcYminus)./(2.*deltaY);	
      
      d2wcdx2 = (wcXplus + wcXminus - 2.*wcX)./(deltaX.^2);
      d2wcdy2 = (wcYplus + wcYminus - 2.*wcX)./(deltaY.^2);
      
      d2wCdx2plus = (wcX2plus - 2.*wcXplus + wcX)./(deltaX.^2);
      d2wCdx2minus= (wcX - 2.*wcXminus + wcX2minus)./(deltaX.^2);      
      d2wCdy2plus = (wcXplusYplus - 2.*wcXplus + wcXplusYminus)./(deltaY.^2);
      d2wCdy2minus= (wcXminusYplus - 2.*wcXminus + wcXminusYminus)./(deltaY.^2);
      
      d3wcdx3 = (d2wCdx2plus - d2wCdx2minus)./(2.*deltaX);      
      d3wcdxdy2=(d2wCdy2plus - d2wCdy2minus)./(2.*deltaX);
      
      %----------------
      %same story for deriving wd(x) 3 times -> UPDATED V2022-03-12
      dwddx = (wdXplus - wdXminus)./(2.*deltaX);
%       dwddy = (wdYplus - wdYminus)./(2.*deltaY);	
      
      d2wddx2 = (wdXplus + wdXminus - 2.*wdX)./(deltaX.^2);
      d2wddy2 = (wdYplus + wdYminus - 2.*wdX)./(deltaY.^2);
      
      d2wDdx2plus = (wdX2plus - 2.*wdXplus + wdX)./(deltaX.^2);
      d2wDdx2minus= (wdX - 2.*wdXminus + wdX2minus)./(deltaX.^2);      
      d2wDdy2plus = (wdXplusYplus - 2.*wdXplus + wdXplusYminus)./(deltaY.^2);
      d2wDdy2minus= (wdXminusYplus - 2.*wdXminus + wdXminusYminus)./(deltaY.^2);
            
      d3wddx3 = (d2wDdx2plus - d2wDdx2minus)./(2.*deltaX);      
      d3wddxdy2=(d2wDdy2plus - d2wDdy2minus)./(2.*deltaX);

	%------------
       %% Prepare equation 2 - bending moment = 0 on the loaded slab (eq. 13.45)
      M(2,:,i) = [d2wadx2 + nu.*d2wady2, d2wbdx2 + nu.*d2wbdy2, 0, 0];
      N(2,:,i) = -d2wdx2 - nu.*d2wdy2;

      %% Prepare equation 3 - bending moment = 0 on the unloaded slab (eq. 13.46)
      M(3,:,i) = [0, 0,d2wcdx2 + nu.*d2wcdy2, d2wddx2 + nu.* d2wddy2];
    %  N(3,:,i) = -d2wdx2C - nu.*d2wdy2C;
        N(3,:,i) = 0;

      %% PREPARE equation 4 - Shear stress equation  
      %%Update: Correct to follow formulation in Deep et al 20 (A, B)
      M(4,1,i) = d3wadx3 + (2-nu).*d3wadxdy2 - (2*g)./(l.^2)*dwadx;
      M(4,2,i) = d3wbdx3 + (2-nu).*d3wbdxdy2 - (2*g)./(l.^2)*dwbdx;
      M(4,3,i) = -1.*(d3wcdx3 + (2-nu).*d3wcdxdy2 - (2*g)/(l.^2).*dwcdx);
      M(4,4,i) = -1.*(d3wddx3 + (2-nu).*d3wddxdy2 - (2*g)/(l.^2).*dwddx);
      N(4,:,i) = -1.*(d3wdx3  + (2-nu).*d3wdxdy2  - (2*g)/(l.^2).*dwdx);

    %     M(4,1,i) = d3wAdx3C + (2-nu).*d3wAdxdy2C - (2*g)/(l^2)*wac;
    %     M(4,2,i) = d3wBdx3C + (2-nu).*d3wBdxdy2C - (2*g)/(l^2)*wbc;
    %     M(4,3,i) = -1*(d3wCdx3C + (2-nu).*d3wCdxdy2C - (2*g)/(l^2)*wcc);
    %     M(4,4,i) = -1*(d3wDdx3C + (2-nu).*d3wDdxdy2C - (2*g)/(l^2)*wdc);
    %     N(4,:,i) = 0;

      %% SOLVE FOR ABCD
      MM = reshape(M(:,:,i),4,4);
      NN = reshape(N(:,:,i),4,1);
      xx = MM\NN;  %% this can be unstable...
%       alternative -> use Gauss-jordan elimination. Slower [10x slower]
%       but less prone to drop alerts.. Doesn't work. Throws error too!
%       xx = gauss_jordan_elim(MM,NN);

      A(i,j) = xx(1);
      B(i,j) = xx(2);
      C(i,j) = xx(3);
      D(i,j) = xx(4);

    end %end for iteration on s variable
end %end for iteration over j
%%that'd be it
end

