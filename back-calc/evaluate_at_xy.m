function fxy = evaluate_at_xy(x,y,XDOMAIN,YDOMAIN,F)
%function fxy = evaluate_at_xy(x,y,XDOMAIN,YDOMAIN,F)
%
%Auxiliary function that evaluates the function F(XDOMAIN,YDOMAIN) at point x,y
%If no exact value for x,y exists in XDOMAIN,YDOMAIN, the nearest value is returned

%release candidate V2022-05-01
%%
fxy = zeros(length(x),length(y));
XDOMAIN = XDOMAIN(:);
YDOMAIN = YDOMAIN(:);

locx = zeros(size(x));
locy = zeros(size(y));
for i = 1:length(x)
    aux = find(abs(x(i)-XDOMAIN) == min(abs(x(i)-XDOMAIN)));
    locx(i) = aux(1);  %stability fix (at times the find command may return two positions - x(i) exactly in between two points)    
end
for j = 1:length(y)
    aux = find(abs(y(j)-YDOMAIN) == min(abs(y(j)-YDOMAIN)));
    locy(j) = aux(1);
end

%fill up the output!

for i = 1:length(x)    
    for j = 1:length(y)
        fxy(i,j) = F(locx(i),locy(j));
    end
end



end
