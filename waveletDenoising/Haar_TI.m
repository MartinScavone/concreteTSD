function [wa,wd] = Haar_TI(x,Levels)

L = length(x);
if nargin<2
    Levels = floor(log2(L));
end
x = x(:);

wa = zeros(L,Levels);
wd = zeros(L,Levels);
wa(:,1) = (x+shift(x,-1))/sqrt(2);
wd(:,1) = (x-shift(x,-1))/sqrt(2);

for i=2:Levels
    wa(:,i) = (wa(:,i-1)+shift(wa(:,i-1),-2^(i-1)))/sqrt(2);
    wd(:,i) = (wa(:,i-1)-shift(wa(:,i-1),-2^(i-1)))/sqrt(2);
end

wa = wa';
wd = wd';