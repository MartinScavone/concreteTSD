function [wa,wb] = compute_wawb(x,y,s,p,b,k,g,l,A,B)
%function [wa,wb] = compute_wawb(x,y,s,p,b,k,g,l,A,B)
%Compute the loaded-finite-slab term of the deflection solution. Following Van
%Cauwelaert, 2004, chap. 15
%
%INPUTS:
    %x,y: coordinates where to evaluate the term [vector]
    %s: wave domain for the fourier transform
    %p: applied pressure (N/m2] 
    %k: modulus of subgrade reaction [N/m3]
    %g: adim parameter related to subgrade's G
    %l: adim parameter related to slab's E and h    
    %A,B: terms for the boundary condition, as solved for each value of s with solveABCD [sizeS x sizeY] 
%
%OUTPUTS:
    %wa, wb(x,y): homogeneous-equation term of the deflection function for
    %the loaded slab
%
%INTERIM SLOW IMPLEMENTATION. ITERATIVELY SOLVE FOR SINGLE LOCATIONS X,Y
%
%release candidate V2022-05-01
%
%%    
wa = zeros(length(x),length(y));
wb = zeros(length(x),length(y));
s = s(:);
%start iterative loop for all points in x,y domain.
for j = 1:length(y)
  w3 = (1./s).*(cos(s.*y(j)./l).*sin(s.*b./l));
  w3 = w3(:);  %force w3 to be a single column like A, B are (Octave makes w3 a row and fails to do element-by-element products afterward, crashing the entire thing!)
  for i = 1:length(x)
%       xi = abs(x(i));
      xi = x(i);
      if g >1
          %case g > 1 - Equations 15.21. Parameters z1, z2(s) are the same as for w(xy) 
          COEF = p /(2*pi*k).*(g.^2-1).^-0.5; 
          z1 = s.^2 + g + sqrt(g^2-1);
          z1 = sqrt(z1);
          z2 = s.^2 + g - sqrt(g^2-1);
          z2 = sqrt(z2);        
          waAux = A(:,j).*w3.*exp(z1.*xi./l);
          wbAux = B(:,j).*w3.*exp(z2.*xi./l);
          %update with Matlab's built-in trapz function to integrate
          wa(i,j) = COEF.*trapz(s,waAux);
          wb(i,j) = COEF.*trapz(s,wbAux);
      elseif g <1 
          %select case g < 1 - Equations 15.24
          COEF = p./(pi.*k).*(1-g^2)^-0.5; 
          alpha = 0.5.*(sqrt((s.^2+g).^2 +1-g^2)+(s.^2 + g));
          alpha = sqrt(alpha);
          beta  = 0.5.*(sqrt((s.^2+g).^2 +1-g^2)-(s.^2 + g));
          beta  = sqrt(beta);        
          waAux = A(:,j).*w3.*(cos(beta./l.*xi)).*exp(alpha.*xi./l);
          wbAux = B(:,j).*w3.*(sin(beta./l.*xi)).*exp(alpha.*xi./l);
          %update with Matlab's built-in trapz function to integrate
          wa(i,j) = COEF.*trapz(s,waAux);
          wb(i,j) = COEF.*trapz(s,wbAux); 
      else
          %case g == 1 - equation 15.23
          COEF = p /(2*pi*k);  
          z = (s.^2 + 1).^0.5;
          waAux = A(:,j).*w3.*exp(z.*xi/l);
          wbAux = B(:,j).*w3.*xi./l.*exp(z.*xi/l);
          %update with Matlab's built-in trapz function to integrate
          wa(i,j) = COEF.*trapz(s,waAux);
          wb(i,j) = COEF.*trapz(s,wbAux);
      end
  end
end

end %endfunction
    
