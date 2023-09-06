function w = compute_w_term(x,y,s,p,a,b,k,g,l)
%function w = compute_w_term(x,y,s,p,a,b,k,g,l)
%
%Compute the inifinte-beam term of the deflection solution. 
%Following Van Cauwelaert, 2004, chap. 15
%
%INPUTS:
    %x,y: coordinates where to evaluate the term [vector]
    %s: wave term (wave domain) used in the Fourier decomposition of the load
    %p: distributed load (Newton/m)
    %a,b: dimensions [length-width] of load [m]
    %k: coefficient of subgrade reaction
    %g: adim parameter related to subgrade's G
    %l: adim parameter related to slab's E and h    
%
%OUTPUTS:
    %w(x,y): infinite-slab term of the deflection function
%compute auxiliary terms alpha and beta - function of wave number s
%
%INTERIM SLOW IMPLEMENTATION. ITERATIVELY SOLVE FOR SINGLE LOCATIONS X,Y
w = zeros(length(x),length(y));
% x1 = find(x<=a);
% x2 = find(x>a);

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
      %w3 = w3(:);
      for i = 1:length(x)
        if x(i) <=a
          %case x<a
          %compute auxiliary terms w1,w2(x,s)
          w1 = z2^-2  .* (2-exp(-z2./l.*(a-x(i)))-exp(-z2./l.*(a+x(i))));
          w2 = z1^-2  .* (2-exp(-z1./l.*(a-x(i)))-exp(-z1./l.*(a+x(i))));
          w1 = w1(:);
          w2 = w2(:);
          %and then calculate the term within the integral
          w_int = (w1-w2).*w3;
        else
          %case x>a
           %compute auxiliary terms w1,w2(x,s)
          w1 = z2^-2  .* (exp(-z2./l.*(x(i)-a))-exp(-z2./l.*(x(i)+a)));
          w2 = z1^-2  .* (exp(-z1./l.*(x(i)-a))-exp(-z1./l.*(x(i)+a)));
          w1 = w1(:);
          w2 = w2(:);
          %and then calculate the term within the integral
          w_int = (w1-w2).*w3;
        end
        %the term I must pass is w_int
        w(i,j) = COEF*(w_int);  
      end  
    end% end double for loop   

elseif g <1 
    %select case g < 1 - Equations 15.15 - 15.16    
    COEF = p /(pi*k);
    %auxiliary coefficients alpha, beta (s)
    alpha = 0.5*(sqrt((s.^2+g).^2 +1-g^2)+(s.^2 + g));
    alpha = sqrt(alpha);
    beta  = 0.5*(sqrt((s.^2+g).^2 +1-g^2)-(s.^2 + g));
    beta  = sqrt(beta);
    
    for j = 1:length(y)
      %the term w3 of the integral only depends on y and s 
      %compute w3 (y,s)
      w3 = cos(s.*y(j)./l).*sin(s.*b./l);
      w3 = w3./(s.*(s.^4+2.*g.*s.^2+1)); 
      %w3 = w3(:);      
      for i = 1:length(x)
        if x(i) <=a
          %case x<a
          %compute auxiliary terms w1,w2(x,s)
          %BUG FOUND V 2022-02-23 -> tHERE'S A PARENTHESIS MISSING IN EXPRESSION FOR W1. It has been parched up in compute_w. m. but
          %not in here.  Dang!
          w1 = (1-g.^2).^-0.5.*exp(-alpha./l.*(a-x(i)));
          w1 = w1.*(sqrt(1-g.^2).*cos((a-x(i)).*beta./l) + (s.^2+g).*(sin((a-x(i)).*beta./l)));
%           w1 = w1(:);
          w2 = (1-g.^2).^-0.5.*exp(-alpha./l.*(a+x(i)));
          w2 = w2.*(sqrt(1-g.^2).*cos((a+x(i)).*beta./l) + (s.^2+g).*(sin((a+x(i)).*beta./l)));
%           w2 = w2(:);
          %and then calculate the term within the integral
          w_int = (2-w1-w2).*w3;
        else
          %case x>a          
          %compute auxiliary terms w1,w2(x,s)
          %BUG FOUND V 2022-02-23 -> tHERE'S A PARENTHESIS MISSING IN EXPRESSION FOR W1. It has been parched up in compute_w. m. but
          %not in here.  Dang!          
          w1 = exp(-alpha./l.*(x(i)-a)).*(1-g^2)^-0.5;
%         w1 = w1.*((1-g.^2).^0.5.*cos((x(i)-a).*beta./l) + (s.^2+g).*sin((x(i)-a).*beta./l));
          w1 = w1.*(sqrt(1-g.^2).*cos((x(i)-a).*beta./l) + (s.^2+g).*(sin((x(i)-a).*beta./l)));
%         w1 = w1(:);
          w2 = exp(-alpha./l.*(x(i)+a)).*(1-g^2)^-0.5;
          w2 = w2.*(sqrt(1-g.^2).*cos((x(i)+a).*beta./l) + (s.^2+g).*(sin((x(i)+a).*beta./l)));
%           w2 = w2(:);
          %and then calculate the term within the integral
          w_int = (w1-w2).*w3;
        end
        %the term I must pass is w_int
        w(i,j) = COEF*(w_int);  
      end  
    end% end double for loop   

else
    %case g == 1 - equation 15.17-18    
    COEF = p /(2*pi*k);
    z = (s.^2 + 1).^0.5;    
    for j = 1:length(y)
      %the term w3 of the integral only depends on y and s 
      %compute w3 (y,s)
      w3 = cos(s.*y(j)./l).*sin(s.*b./l);
      w3 = w3./(s.*((s.^2+1)^2));
      %w3 = w3(:);    
      for i = 1:length(x)
        if x(i) <=a
          %case x<a
          %compute auxiliary terms w1,w2(x,s)
          w1 = exp(-z./l.*(a-x(i)));
          w1 = w1.*(2 + sqrt(1+s.^2).*(a-x(i))./l);
          w2 = exp(-z./l.*(a+x(i)));
          w2 = w2.*(2 + sqrt(1+s.^2).*(a+x(i))./l);
          w1 = w1(:);
          w2 = w2(:);
          %and then calculate the term within the integral
          w_int = (4-w1-w2).*w3;
        else
          %case x>a          
          %compute auxiliary terms w1,w2(x,s)
          w1 = exp(-z./l.*(x(i)-a));
          w1 = w1.*(2 + sqrt(1+s.^2).*(x(i)-a)./l);
          w2 = exp(-z./l.*(x(i)+a));
          w2 = w2.*(2 + sqrt(1+s.^2).*(x(i)+a)./l);
          w1 = w1(:);
          w2 = w2(:);
          %and then calculate the term within the integral
          w_int = (w1-w2).*w3;
        end
        %the term I must pass is w_int
        w(i,j) = COEF*(w_int);  
      end  
    end% end double for loop     
end

end %endfunction    