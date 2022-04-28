function y = softShrinking(x,lambda)
% function y = softShrinking(x,t)
% This function will apply the 'soft shrinking rule' to the vector x,
% guided by the value of 't'
%See the paper by Katicha et al. (2013) on wavelet denoising for info on
%what is this about.
%
%Input: x: vector of data
%       t: positive scalar (not necessary integer)
%Output:y: vector the size of x, filtered as per the soft shrinking rule.
%
%Release candidate v2022-05-01

%% integrity check. t must be a positive scalar

integrityCheck = lambda>=0;
if ~integrityCheck
   warning('softShrinking:: parameter "t" is smaller than 0')
   return
end

%Prep work: Compatibility for the both BPD and Rew-L1-decompos
%BP can work with a scalar lambda, whereas Reweighed L1 decomposition needs
%lambda as a vector (cause each entry of x may have it's own lambda(x) =
%lambda*w(x)

%check that if lambda enters as a scalar it gets converted to a vector the
%size of x (so that it doesn't crash below when doing the soft-shrinking)

if length(lambda) == 1
    lambda = lambda.*ones(size(x));
else
    %do nothing
end


%% apply the rule - treat separately those values of x with absolute value greater or less than t
y = zeros(size(x));

zone1 = find(abs(x)<=lambda);
zone2 = find(abs(x)>lambda);

if ~isempty(zone1)
    %the positions of y told by zone1 must be filled up with zeros
    y(zone1) = zeros(size(zone1));
end
if ~isempty(zone2)
    %the positions of y told by zone2 must be filled up with y =
    %x-t*signum(x). Note, sign(x) = signum function for x
    y(zone2) = x(zone2) - sign(x(zone2)).*lambda(zone2);
end

