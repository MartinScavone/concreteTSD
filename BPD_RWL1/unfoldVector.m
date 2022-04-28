function w = unfoldVector(vector,n)
%function w = unfoldVector(vector,n)
%function to unfold a vector 'vector' of any given length to a length of n
%by adding extra elements at both ends in a symmetric way (mirroring the
%elements of 'vector' over the two extremes).
%NOTE: If I need to add an even number of elements, it will add an equal size of entries both ahead and after the vector
%Otherwise, it will add one more entry to the beginning of the vector
%
%e.g.
%If vector is [1 2 3 4 5 6] and n = 10, then w = [3 2 1 2 3 4 5 6 5 4]
%If vector is [1 2 3 4 5]   and n = 10, then w = [4 3 2 1 2 3 4 5 4 3]
%
%Release candidate v2022-05-01
%
%V0.1 2019-08-16
%Update: 1) Don't repeat the extreme values of 'vector'.
%        2) allow for n to also be an odd number.  (remove the "isItEven" restriction for n)

%% 1) calculate how many values do I need to add.

len = length(vector);
[a,b] = size(vector);

ext = n-len;
% check if ext is even or odd and so define the number of elements to add
% to the beginning and end of 'vector' to assemble w [called nw1 and nw2
% respectively]

isItEven = ext/2 == floor(ext/2);
if isItEven
   %add an equal number to nw1 and nw2
    nw1 = floor(ext/2);
    nw2 = floor(ext/2);   
else
    %add one more element to nw1 than nw2
    nw1 = floor(ext/2)+1;
    nw2 = floor(ext/2);
end

%% 2) Now build the w vector. Assemble w1, w2 as row vectors. Transpose at the end to match the shape of 'vector' if necessary

w1 = vector(nw1+1:-1:2);
w2 = vector(end-1:-1:end-nw2);

%Finally adapt the shape of the output vector 'w' to that of 'vector'
if a>=b
    %vector is a column -> transpose w1 and w2 and merge with vector
    w = [w1; vector; w2];
else
    %vector is also a row -> paste all pieces together
    w = [w1 vector w2];
end

   
end  %% endfunction
