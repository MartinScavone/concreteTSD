function out = Wth(x,p,th)

if nargin<3
    th = sqrt(2*log(length(x)));
    if nargin<2
        p = 2;
    end
end

% p = max(p,1);

if p>1e4
    out = x;
    out(abs(out)<=th) = 0;
else
    s = sign(x);
    out = max(abs(x).*(1-(th./abs(x)).^p),0);
    out = s.*out;
end