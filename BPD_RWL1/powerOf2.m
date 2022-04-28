function [m,expo] = powerOf2(x)
%function m = powerOf2(x)
%This function is to locate the closest power of 2 that is greater than x
%OUTPUT:
% m = smallest power of 2 number that is greater than x
% expo= exponent to raise 2 to get m  (expo = log2(m))
%
%Release candidate v2022-05-01

%Dummy coding, it will work by scanning all powers of 2 until it surpasses
%the value of X
%If by any chance it doesn't find anything, it will return m = -1 and expo
%= 0
m = -1;
expo = 0;
while m == -1
    %do the power of 2
    tempM = 2^expo;
    %try if the power of 2 (tempM) is greater than input number x
    if tempM >= x
        m = tempM;
    else
        %if not, try the next power of 2
        expo = expo+1;
    end
    
end