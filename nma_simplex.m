function tab = nma_simplex(A,b,c,debug)
%function [A,b,c]=nma_simplex(A,b,c)
%This function implments the simplex matrix algorithm.
%It accepts A_eq and b_eq and c as defined in standard 
%documentation and generates all the simplex tableaus, and 
%returns the final tableau which the user can read from it the 
%minimum value of the objective function and the optimal x vector 
%directly.
%
%It runs both phase one and phase two automatically.
%
%The input is
%
%A: This is the Ax=b matrix. This is for simplex standard 
%   form only. The caller must convert all inequalites to 
%   equalities first by using slack and suprluse variables. This 
%   is what is called the Aeq matrix in Matlab documenation. 
%   This function does not support Ax<b form. A has to be in 
%   standard form
%
%b: Vector. This is the right hand side of Ax=b.
%
%c: Vector. This is from minimize  J(x) = c'x. As defined in 
%   standard Matlab documentations.
%
%debug: flag. Set to true to see lots of internal steps.
%
%Returns:
%
%This function returns the final tableau. It has the form
%
%  [ A | b  ]
%  [ c | J  ]
%
% Version 5/12/2016
% by Nasser M. Abbasi
% Free for use.

validate_input(A,b,c);

[A,b] = make_phase_one(A,b,debug);
tab   = simplex(A,b,c,debug,'phase two');
end
%==========================
function [A,b] = make_phase_one(A,b,debug)
[m,n]              = size(A);
tab                = zeros(m+1,n+m+1);
tab(1:m,1:n)       = A;
tab(end,n+1:end-1) = 1;
tab(1:m,end)       = b(:);
tab(1:m,n+1:n+m)   = eye(m);

if debug
    fprintf('>>>>Current tableau [phase one]\n');
    disp(tab);
end

for i = 1:m %now make all entries in bottom row zero
    tab(end,:) = tab(end,:)-tab(i,:);
end

tab = simplex(tab(1:m,1:n+m),tab(1:m,end),tab(end,1:n+m),...
                                             debug,'phase one');
%if tab(end,end) ~=0
%   error('artificial J(x) is not zero at end of phase one.');
%end

A = tab(1:m,1:n);
b = tab(1:m,end);

end
%=================================
function tab = simplex(A,b,c,debug,phase_name)
[m,n]        = size(A);
tab          = zeros(m+1,n+1);
tab(1:m,1:n) = A;
tab(m+1,1:n) = c(:);
tab(1:m,end) = b(:);

keep_running = true;
while keep_running
    if debug
        fprintf('***********************\n');
        fprintf('Current tableau [%s] \n',phase_name);
        disp(tab);
    end
    
    if any(tab(end,1:n)<0)%check if there is negative cost coeff.
        [~,J] = min(tab(end,1:n)); %yes, find the most negative
        % now check if corresponding column is unbounded
        if all(tab(1:m,J)<=0) 
          error('problem unbounded. All entries <= 0 in column %d',J);
        %do row operations to make all entries in the column 0 
        %except pivot
        else 
            pivot_row = 0;
            min_found = inf;
            for i = 1:m
                if tab(i,J)>0
                    tmp = tab(i,end)/tab(i,J);
                    if tmp < min_found
                        min_found = tmp;
                        pivot_row = i;
                    end
                end
            end
            if debug
                fprintf('pivot row is %d\n',pivot_row);
            end
            %normalize
            tab(pivot_row,:) = tab(pivot_row,:)/tab(pivot_row,J); 
            %now make all entries in J column zero.
            for i=1:m+1
                if i ~= pivot_row
                    tab(i,:)=tab(i,:)-sign(tab(i,J))*...
                                 abs(tab(i,J))*tab(pivot_row,:);
                end
            end
        end
        if debug  %print current basic feasible solution
            fprintf('current basic feasible solution is\n');
            disp(get_current_x());
        end
    else
        keep_running=false;       
    end
end

    %internal function, finds current basis vector
    function current_x = get_current_x()
        current_x = zeros(n,1);
        for j=1:n
            if length(find(tab(:,j)==0))==m
                idx= tab(:,j)==1;
                current_x(j)=tab(idx,end);
            end
        end
    end
end
%==========================
function validate_input(A,b,c)
if ~ismatrix(A)
    error('A must be matrix');
end

if ~isvector(b)
    error('b must be vector');
end
if ~isvector(c)
    error('c must be vector');
end

[m,n]=size(A);
if rank(A) <m
    error('Rank A must be equal to number of rows in A');
end

if length(b) ~= m
    error('b must have same size as number of rows of A');
end
if length(c) ~= n
    error('c must have same size as number of columns of A');
end
end