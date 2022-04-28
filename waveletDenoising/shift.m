function y = shift(x,shift_size)
%
% y = shift(x,shift_size) cicularly shifts the rows in matrix x by
% shif_size positions. If x is a vector, the shift is performed on the
% elements of x

[r,c] = size(x);
id = 1:length(x);

if shift_size>0
    id = [id(end-shift_size+1:end) id(1:end-shift_size)];
elseif shift_size<=0
    id = [id(1-shift_size:end) id(1:-shift_size)];
end

if r==1 || c==1
    y = x(id);
else
    y = x(id,:);
end