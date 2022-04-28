function vector = foldBackVector(w,len)
%function vector = foldBackVector(w,len)
%This function rolls back the effect of 'unfoldVector': a vector of length
%'len' is recovered from the unfolded/mirrored vector w.
%
%See the info for 'unfoldVector' for details on how 'w' is constructed from
%'vector'
%
%Release candidate v2022-05-01

%initialize output
[a,b] = size(w);
if a>=b
    %w is a column vector
    vector = zeros(len,1);
else
    vector = zeros(1,len);
end

%UPDATE v2019-10-27: I had sth. wrong with the location of the first and
%last elements of the soruce vector "vector"
%re'do the entire thing.
addedElements = length(w)-len;  %all the elements that were added when unfolding.
addedToTheEnd = floor(addedElements/2);   %the unfolding adds one less element at the end if the source vector is odd-size, or half and half if it's even.
addedBeginnin = addedElements - addedToTheEnd;

%now extract vector from w
vector = w(addedBeginnin+1:addedBeginnin+len);
