function [wc,wd] = compute_wcwd_term(x,y,s,p,b,k,g,l)
%function [wc,wd] = compute_wcwd_term(x,y,s,p,b,k,g,l)
%
%Axuliary function for the solveABCD function. Compute the term of the WcWd terms
%for a single value of the wave number s.
%Following Van %Cauwelaert, 2004, chap. 15
%
%INPUTS:
    %x,y: coordinates where to evaluate the term [vector]    
    %s: wave domain for the fourier transform
    %p: applied pressure (N/m2]
    %k: modulus of subgrade reaction [N/m3]
    %g: adim parameter related to subgrade's G
    %l: adim parameter related to slab's E and h  
%
%OUTPUTS:
    %wc, wd(x,y): homogeneous-equation term of the deflection function for
    %the unloaded slab
%
%release candidate V2022-05-01
    
%%
wc = zeros(length(x),length(y));
wd = zeros(length(x),length(y));
%%

for j = 1:length(y)
  w3 = (1./s).*(cos(s.*y(j)./l).*sin(s.*b./l));
  w3 = w3(:);
  for i = 1:length(x)
    %%start case for wa(x,y), wb (x,y)
    if g >1
        %case g > 1 - Equations 15.21. Parameters z1, z2(s) are the same as for w(xy)
        COEF = p /(2*pi*k).*(g.^2-1).^-0.5;         
        z1 = s.^2 + g + sqrt(g^2-1);
        z1 = sqrt(z1);
        z2 = s.^2 + g - sqrt(g^2-1);
        z2 = sqrt(z2);
        
        wc(i,j) = COEF.*w3.*exp(-z1.*x(i)./l);
        wd(i,j) = COEF.*w3.*exp(-z2.*x(i)./l);
                       
    elseif g <1 
        COEF = p /(pi*k)*(1-g^2)^-0.5;
        %auxiliary coefficients alpha, beta (s)
        alpha = 0.5*(sqrt((s.^2+g).^2 +1-g^2)+(s.^2 + g));
        alpha = sqrt(alpha);
        beta = 0.5*(sqrt((s.^2+g).^2 +1-g^2)-(s.^2 + g));
        beta = sqrt(beta);
        wc(i,j) = COEF.*w3.*(cos(beta/l.*x(i))).*exp(-alpha.*x(i)/l);
        wd(i,j) = COEF.*w3.*(sin(beta/l.*x(i))).*exp(-alpha.*x(i)/l);        
        
    else
        COEF = p /(2*pi*k);
        z = (s.^2 + 1).^0.5;
        wc(i,j) = COEF.*w3.*exp(-z.*x(i)/l);
        wd(i,j) = COEF.*w3.*x(i)./l.*exp(-z.*x(i)/l);     
        
    end
  end
end

end %endfunction    

