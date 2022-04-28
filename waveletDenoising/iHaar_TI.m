function x = iHaar_TI(wa,wd)

Levels = size(wd,1);
x = (wa(end,:)+wd(end,:)+shift((wa(end,:)-wd(end,:))',2^(Levels-1))')/2/2*sqrt(2);

for i=Levels-1:-1:1
    x = (x+wd(i,:)+shift((x-wd(i,:))',2^(i-1))')/2/2*sqrt(2);
end