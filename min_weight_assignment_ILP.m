function [i_max,weight_1,fval] = min_weight_assignment_ILP(path_length,sharedpath_length)
% compute the link weight assignment consistent with given path_length and
% sharedpath_length that has the minimum total weight
% output: weight(i) is the weight of link (sum weight of all links)
% traversed exactly by paths in Acal{i}.

% Note: path_length must be in the order of nchoosek(1:n,1);
% sharedpath_length must be in the order of nchoosek(1:n,2).
    
n = length(path_length); % #paths
for i_max=1:1:2^n-1
% use: x = linprog(f,A,b,Aeq,beq,lb,ub)
    combos = combntns(1:(2^n-1),i_max);
    f = ones(2^n-1,1);
    A=[];
    b=[];
    Aeq_f=zeros(n*(n-1)/2+n,2^n-1); % B(i,j)=1 iff Acal(i) subseteq Acal(j)
    % compute the coefficient matrix:
    Acal = construct_A(n);
    for i=1:n*(n-1)/2+n
        for j=1:2^n-1
            Aeq_f(i,j) = all(ismember(Acal{i}, Acal{j}));
        end
    end
    beq_f = zeros(n*(n-1)/2+n,1);
    beq_f(1:n) = path_length;
    beq_f(n+1:end) = sharedpath_length;
    lb = zeros(2^n-1,1);
    ub = [];
    
    for j_loop=1:1:nchoosek(2^n-1,i_max)
        Aeq_add=zeros(i_max,2^n-1);
        for k=1:1:i_max
           Aeq_add(k,combos(j_loop,k))=1; 
        end
        beq_add=zeros(i_max,1);
        Aeq = [Aeq_f; Aeq_add];
        beq = [beq_f; beq_add];
        [weight,fval,exitflag] = linprog(f,A,b,Aeq,beq,lb,ub);
        if(exitflag==1) 
            i_max=i_max+1;
            weight_1=weight;
            break;
        end
        
        if(j_loop==nchoosek(2^n-1,i_max))
            i_max=i_max-1;
            return;
        end
    end
end
if(i_max==2^n-1) 
    return;
end

end

