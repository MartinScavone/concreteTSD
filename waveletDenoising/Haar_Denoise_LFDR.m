function [yd,swa,swd] = Haar_Denoise_LFDR(yn,fdr,s)

if nargin<3
    s = 1.4826*median(abs(diff(yn)-median(diff(yn))))/sqrt(2);
    if nargin<2
        fdr = 0.5;
    end
end
L = length(yn);

[swa,swd] = Haar_TI(yn);

for i=1:size(swd,1)
    p = 2*(1-normcdf(abs(swd(i,:)),0,s));
    pth = vfdr(p(:),fdr);
    th = max(abs(swd(i,p>=pth)));
    if isempty(th)
        th = max(abs(swd(i,:)))+1;
    end
%     swd(i,:) = scad_th(swd(i,:),th,sqrt(2*log(L))*s);
    swd(i,:) = Wth(swd(i,:),1000000,th);
end

yd = iHaar_TI(swa,swd);
yd = yd(:);