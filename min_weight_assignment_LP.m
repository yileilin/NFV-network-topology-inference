function [weight,fval,exitflag] = min_weight_assignment_LP(path_length,sharedpath_length)
% compute the link weight assignment consistent with given path_length and
% sharedpath_length that has the minimum total weight
% output: weight(i) is the weight of link (sum weight of all links)
% traversed exactly by paths in Acal{i}.

% Note: path_length must be in the order of nchoosek(1:n,1);
% sharedpath_length must be in the order of nchoosek(1:n,2).
n = length(path_length); % #paths
% use: x = linprog(f,A,b,Aeq,beq,lb,ub)
f = ones(2^n-1,1);
A=[];
b=[];
Aeq=zeros(n*(n-1)/2+n,2^n-1); % B(i,j)=1 iff Acal(i) subseteq Acal(j)
% compute the coefficient matrix:
Acal = construct_A(n);
for i=1:n*(n-1)/2+n
    for j=1:2^n-1
        Aeq(i,j) = all(ismember(Acal{i}, Acal{j}));
    end
end
beq = zeros(n*(n-1)/2+n,1);
beq(1:n) = path_length;
beq(n+1:end) = sharedpath_length;
lb = zeros(2^n-1,1);
ub = [];
%x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub)

temp=nma_simplex(Aeq,beq,f,false);
weight=zeros(2^n-1,1);
for i=1:1:2^n-1
    sum_1=0;
    sum_0=0;
    index_1=0;
    for j=1:1:n*(n-1)/2+n
        if temp(j,i)==1
            sum_1=sum_1+1;
            index_1=j;
        end
        if temp(j,i)==0
            sum_0=sum_0+1;
        end
    end
    if sum_1==1 && sum_0==n*(n-1)/2+n-1
        weight(i,1)=temp(index_1,2^n);
    end
end

fval=f'*weight;
exitflag=1;
%[weight,fval,exitflag] = linprog(f,A,b,Aeq,beq,lb,ub);

end


