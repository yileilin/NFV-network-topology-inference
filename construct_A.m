function [Acal] = construct_A(n)
% construct an array of sets Acal in the following order: all the size-1
% subsets of {1,...,n}, all the size-2 subsets of {1,...,n},... set
% {1,...,n} itself.
Acal = cell(1,2^n-1);
j=0;
for i=1:n
    Atmp = nchoosek(1:n,i);
    for k=1:nchoosek(n,i)
        Acal{j+k} = Atmp(k,:);
    end
    j = j + nchoosek(n,i);
end
end