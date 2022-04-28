function w = compute_w(x,y,s,p,a,b,k,g,l)
%function w = compute_w(x,y,s,p,a,b,k,g,l)
%Compute the infinite-beam term of the deflection solution. Following Van
%Cauwelaert, 2004, chap. 15
%
%INPUTS:
    %x,y: coordinates where to evaluate the term [vector]
    %s: wave term (wave domain) used in the Fourier decomposition of the load
    %p: pressure of the applied load (Newton/m2)
    %a,b: dimensions [half length - half width] of load [m]
    %k: coefficient of subgrade reaction [N/m3]
    %g: adim parameter related to subgrade's G
    %l: adim parameter related to slab's E and h    
%
%OUTPUTS:
    %w(x,y): infinite-slab term of the deflection function
%INTERIM SLOW IMPLEMENTATION. ITERATIVELY SOLVE FOR SINGLE LOCATIONS X(i),Y(j)
%
%
%release candidate V2022-05-01

%%
w = zeros(length(x),length(y));
%stability check -in case the wave domain has a 0, remove it (prevent a division by zero)
if s(1) ==0
  s(1) = eps;
end
%% 
if g >1  
    %case g > 1 - Equations 15.13 and 15.14
    COEF = p /(2*pi*k).*(g^2-1).^-0.5;
    %auxiliary coefficients z1, z2 (s)
    z1 = s.^2 + g + sqrt(g^2-1);
    z1 = sqrt(z1);
    z2 = s.^2 + g - sqrt(g^2-1);
    z2 = sqrt(z2);
    for j = 1:length(y)
      %the term w3 of the integral only depends on y and s
      w3 = cos(s.*y(j)./l).*sin(s.*b./l);
      w3 = w3./s;
      w3 = w3(:);
      for i = 1:length(x)
          xi = abs(x(i));
          if xi <=a
              %case abs(x)<a [within the wheel print]
              %compute auxiliary terms w1,w2(x,s)
              w1 = (z2.^-2) .* (2-exp(-z2./l.*(a-xi))-exp(-z2./l.*(a+xi)));
              w2 = (z1.^-2) .* (2-exp(-z1./l.*(a-xi))-exp(-z1./l.*(a+xi)));
              %and then calculate the term within the integral
              w_int = (w1-w2).*w3;
          else
              %case abs(x)>a [outside the wheel print]
              %compute auxiliary terms w1,w2(x,s)
              w1 = (z2.^-2) .* (exp(-z2./l.*(xi-a))-exp(-z2./l.*(xi+a)));
              w2 = (z1.^-2) .* (exp(-z1./l.*(xi-a))-exp(-z1./l.*(xi+a)));
              %and then calculate the term within the integral
              w_int = (w1-w2).*w3;
          end
          %now i must integrate the w_int over S - sum all terms!
%           w(i,j) = COEF.*sum(w_int).*deltaS; 
%           <--improved using the built-in numerical integration function "trapz"
          w(i,j) = (COEF).*trapz(s,w_int);
      end  
    end% end double for loop   

elseif g <1 
    %select case g < 1 - Equations 15.15 - 15.16    
    COEF = p /(pi*k);
    %auxiliary coefficients alpha, beta (s)
    alpha = 0.5.*(sqrt((s.^2+g).^2 +1-g^2)+(s.^2 + g));
    alpha = sqrt(alpha);
    beta  = 0.5.*(sqrt((s.^2+g).^2 +1-g^2)-(s.^2 + g));
    beta  = sqrt(beta);
    
    for j = 1:length(y)
      %the term w3 of the integral only depends on y and s 
      %compute w3 (y(j),s)
      w3 = cos(s.*y(j)./l).*sin(s.*b./l);
      w3 = w3./(s.*(s.^4+2*g.*s.^2+1));
%       w3 = w3(:);      
      for i = 1:length(x)
          xi = abs(x(i));
          if xi <=a
              %case abs(x)<a
              %compute auxiliary terms w1,w2(x,s)
              w1 = (1-g.^2).^-0.5.*exp(-1.*alpha./l.*(a-xi));
              w1 = w1.*((sqrt(1-g.^2)).*(cos((a-xi).*beta./l)) + (s.^2+g).*(sin((a-xi).*beta./l)));
              w2 = (1-g.^2).^-0.5.*exp(-1.*alpha./l.*(a+xi));
              w2 = w2.*((sqrt(1-g.^2)).*(cos((a+xi).*beta./l)) + (s.^2+g).*(sin((a+xi).*beta./l)));
              %and then calculate the term within the integral
              w_int = (2-w1-w2).*w3;
          else
              %case abs(x)>a
              %compute auxiliary terms w1,w2(x,s)
              w1 = exp(-1.*alpha./l.*(xi-a));
              w1 = w1.*((sqrt(1-g.^2)).*(cos((xi-a).*beta./l)) + (s.^2+g).*(sin((xi-a).*beta./l)));
              w2 = exp(-1.*alpha./l.*(xi+a));
              w2 = w2.*((sqrt(1-g.^2)).*(cos((xi+a).*beta./l)) + (s.^2+g).*(sin((xi+a).*beta./l)));
              %and then calculate the term within the integral
              w_int = (w1-w2).*w3./sqrt(1-g^2);
          end
          %now i must integrate the w_int over S - sum all terms!
          %w(i,j) = (COEF).*sum(w_int).*deltaS;
          %<--improved using the built-in numerical integration function "trapz"
          w(i,j) = (COEF).*trapz(s,w_int);   
      end
    end% end double for loop

else
    %case g == 1 - equation 15.17-18    
    COEF = p /(2*pi*k);
    z = (s.^2 + 1).^0.5;    
    for j = 1:length(y)
      %the term w3 of the integral only depends on y and s 
      %compute w3 (y,´s)
      w3 = cos(s.*y(j)./l).*sin(s.*b./l);
      w3 = w3./(s.*((s.^2+1)^2));    
      w3 = w3(:);
      for i = 1:length(x)
          xi = abs(x(i));
          if xi <=a
              %case x<a
              %compute auxiliary terms w1,w2(x,s)
              w1 = exp(-z./l.*(a-xi));
              w1 = w1.*(2 + sqrt(1+s.^2).*(a-xi)./l);
              w2 = exp(-z./l.*(a+xi));
              w2 = w2.*(2 + sqrt(1+s.^2).*(a+xi)./l);
              %and then calculate the term within the integral
              w_int = (4-w1-w2).*w3;
          else
              %case x>a
              %compute auxiliary terms w1,w2(x,s)
              w1 = exp(-z./l.*(xi-a));
              w1 = w1.*(2 + sqrt(1+s.^2).*(xi-a)./l);
              w2 = exp(-z./l.*(xi+a));
              w2 = w2.*(2 + sqrt(1+s.^2).*(xi+a)./l);
              %and then calculate the term within the integral
              w_int = (w1-w2).*w3;
          end
          %now i must integrate the w_int over S - sum all terms!
%           w(i,j) = COEF.*sum(w_int).*deltaS;
%           <--improved using the built-in numerical integration function "trapz"
          w(i,j) = (COEF).*trapz(s,w_int);
      end
    end% end double for loop
end

end %endfunction    