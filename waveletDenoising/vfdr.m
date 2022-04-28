function [pth,q,a,t1t2error,qa] = vfdr(pval,locfdr,a)
%
% function [pth,q] = vfdr(pval,a)
% variable false discovery rate (FDR) procedure to minimize classification
% error
% inputs:
% - pval: the p values of th eobservations
% - a (optional): proportion of observations that come from the null
%                 distribution. If a is not provided, it is estimated from
%                 the data
% Outputs:
% - pth: the p value threshold that results in minimixzing the
%        classification error; p values lower than pth are estimated as not
%        coming from the null distribution while p values larger than pth
%        are estimated as coming from the null distribution. p values that
%        are equal to pth are are estimated as not coming from the null
%        distribution is their resulting qval is less than 0.5, otherwise
%        they are estimated as coming from the null distribution
% - q: the q value corresponding to pth. If q<0.5 then pth is estimated as
%      not coming from the null distribution
% - a: estimated proportion of measurements that come from the null
%      distribution
% - t1t2error: the estimated total classification error for every possible
%              threshold

m = length(pval);
p = sort(pval);
k = 1:m; k = k(:);
if nargin<3
    q = p.*m./k; %q = min(q,1);
    qa = p.*(m+1-k)./(k-k.*p); qa = min(qa,1);
%     a = sum(p<0.75)/m/0.75;
%     a = min(q(round(end/2):end)\qa(round(end/2):end),1);
    a = min(sum(p>=0.5)/(0.5*length(p)),1);
%     a = min(max(a,0.0001),1);
    if nargin<2
        locfdr = 0.5;
    end
end

qval = a*p.*m./k;

for j=length(qval)-1:-1:1
    qval(j) = min(qval(j),qval(j+1));
end

% t1t2error = -k.*(1-2*qval);
t1t2error = -(locfdr*k/m-a*p);
[~,ix] = min(t1t2error);

pth = p(ix);
q = qval(ix);

count = sum(p<=pth);
if count==1
    if q>=0.5
        pth = -1;
    end
end