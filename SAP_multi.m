%TODO
%dest_table  node-i(i=1:1:n_dest)
%cate_table  
function [Aeq,Beq,Aneq,Bneq,result,exit]=SAP_multi(weight,path_table,Aneq,Bneq,i_neq_left)
% n_category=n_path*(n_path-1)/2+n_path;
% l_max=N_VNF+3*n_path*n_path;
% n_dummy=2*n_path*n_path;
% f_total=1+n_dest+N_VNF+n_dummy; %s,d,chain,dummy
% 
%  %alphabet
% x(n_path,l_max,f_total);
% delta_k(1,n_dummy);
% delta_ij(f_total,f_total);
% g1(f_total,f_total,l_max,l_max,n_path);
% g2(f_total,f_total,n_path);
% g3(f_total,f_total,n_category);

%% main

%define some variables
global n_path
global l_max
global n_dest
global f_total
global n_dummy
global n_category
global n_source

N_VNF=5;

% n_path=6;
% l_max=12;
% n_dummy=5;
% n_dest=3;
% n_source=2;

% n_path=4;
% l_max=11;
% n_dummy=5;
% n_dest=2;
% n_source=2;

n_path=2;
l_max=9;
n_dummy=2;
n_dest=1;
n_source=2;

f_total=n_source+n_dest+N_VNF+n_dummy;

%dest_aug
for i=1:n_source
    dest_aug(((i-1)*n_dest+1):(i*n_dest))=(n_source+1):(n_source+n_dest);
end

path_table_new=cell(1,n_path);
for j=1:n_source
    for i=1:n_dest
        path_table_new{1,(j-1)*n_source+i}=path_table{j,i};
    end
end
path_table=path_table_new;
%path_table_aug
path_table_aug=cell(1,n_path);
for k=1:n_source
    for i=1:1:n_dest
        path_table_aug{(k-1)*n_dest+i}=[];
        path_table_aug{(k-1)*n_dest+i}=[path_table_aug{(k-1)*n_dest+i},k];
        for j=1:1:length(path_table{(k-1)*n_dest+i})
            path_table_aug{(k-1)*n_dest+i}=[path_table_aug{(k-1)*n_dest+i},path_table{(k-1)*n_dest+i}(j)+2+n_dest];
        end
        path_table_aug{(k-1)*n_dest+i}=[path_table_aug{(k-1)*n_dest+i},dest_aug((k-1)*n_dest+i)];
    end
end

%A
B=construct_A(n_path);
A=cell(1,1);
index=1;
for i=1:1:2^n_path-1
    if weight(i)>1e-8
        A{index}=B{i};
        index=index+1;
    end
end
n_category=length(A);

total_v=n_path*l_max*f_total+n_dummy+f_total*f_total+ ...,
        f_total*f_total*l_max*l_max*n_path+f_total*f_total*n_path+ ...,
        f_total*f_total*n_category;

Aeq=sparse(1,total_v);
Beq=sparse(1,1);

i_eq_left=1;
i_eq_right=1;

%(8)
% for i=1:1:n_path
%     for j=1:1:l_max
%         Aneq(i_neq_left,x(i,j,1):x(i,j,f_total))=1;
%         i_neq_left=i_neq_left+1;
%     end
% end
% 
% Bneq(1:(n_path*l_max),1)=1;
% i_neq_right=n_path*l_max+1;

%(11)
for i=1:1:n_path
    for f=1:1:n_source+n_dest+N_VNF
        Aeq(i_eq_left,x(i,1:l_max,f))=1;
        i_eq_left=i_eq_left+1;  
    end
end
for i=1:1:n_path
%     if (i<3)
%         Beq(i_eq_right,1)=1;
%         i_eq_right=i_eq_right+1;
%         Beq(i_eq_right,1)=0;
%         i_eq_right=i_eq_right+1;
%     else
%         Beq(i_eq_right,1)=0;
%         i_eq_right=i_eq_right+1;
%         Beq(i_eq_right,1)=1;
%         i_eq_right=i_eq_right+1;        
%     end
%     
%     for f=3:1:2+n_dest
%         if(ismember(f,dest_aug))
%             if i==f-1
%                 Beq(i_eq_right,1)=1;
%             end
%             i_eq_right=i_eq_right+1;
%         end
%     end
    for f=1:1:n_source+n_dest+N_VNF
        if ismember(f,path_table_aug{1,i}) %VNF in path i
            Beq(i_eq_right,1)=1;
        else
            Beq(i_eq_right,1)=0;
        end
        i_eq_right=i_eq_right+1;
    end
end
        
% %(13)
% for f=2+n_dest+N_VNF:1:f_total
%     for i=1:1:n_path
%         for j=1:1:l_max
%             Aneq(i_neq_left,x(i,j,f))=1;
%             Aneq(i_neq_left,delta_k(f-1-n_dest-N_VNF))=-1;
%             i_neq_left=i_neq_left+1;
%         end
%     end
% end
% 
% Bneq(i_neq_right:(i_neq_left-1),1)=0;
% i_neq_right=i_neq_right+1;

% %(15)
% for i=1:1:n_path
%     for f1=1:1:f_total
%         for f2=1:1:f_total
%             for j1=1:1:l_max
%                 for j2=1:1:l_max
%                     Aneq(i_neq_left,g1(f1,f2,j1,j2,i))=1;
%                     Aneq(i_neq_left,x(i,j1,f1))=-1;
%                     i_neq_left=i_neq_left+1;
%                     Aneq(i_neq_left,g1(f1,f2,j1,j2,i))=1;
%                     Aneq(i_neq_left,x(i,j2,f2))=-1;
%                     i_neq_left=i_neq_left+1;
%                     Aneq(i_neq_left,g1(f1,f2,j1,j2,i))=-1;
%                     Aneq(i_neq_left,x(i,j1,f1))=1;
%                     Aneq(i_neq_left,x(i,j2,f2))=1;
%                     Bneq(i_neq_left,1)=1;
%                     i_neq_left=i_neq_left+1;             
%                 end
%             end
%         end
%     end
% end
% i_neq_right=i_neq_left;

%(18)
for i=1:1:n_path
    for p=1:1:length(path_table_aug{i})-1
        f1=path_table_aug{i}(p);
        f2=path_table_aug{i}(p+1);
        for j1=1:1:l_max-1
            for j2=j1+1:1:l_max
                Aeq(i_eq_left,g1(f1,f2,j1,j2,i))=1;               
            end
        end
        i_eq_left=i_eq_left+1;
    end
end
Beq(i_eq_right:(i_eq_left-1),1)=1;
i_eq_right=i_eq_left;

% %(16)
% for i=1:1:n_path
%     for f1=1:1:f_total
%         for f2=1:1:f_total
%             for j=1:1:l_max-1
%                 Aneq(i_neq_left,g1(f1,f2,j,j+1,i))=1;
%                 Aneq(i_neq_left+1,g1(f1,f2,j,j+1,i))=-1;
%             end
%             Aneq(i_neq_left,g2(f1,f2,i))=-l_max;
%             Aneq(i_neq_left+1,g2(f1,f2,i))=1;
%             i_neq_left=i_neq_left+2;
%         end
%     end
% end
% Bneq(i_neq_right:(i_neq_left-1),1)=0;
% i_neq_right=i_neq_left;

% %(20)
% for i=1:1:n_path
%     for f1=1:1:f_total
%         for f2=1:1:f_total
%             Aneq(i_neq_left,g2(f1,f2,i))=1;
%             Aneq(i_neq_left,delta_ij(f1,f2))=-1;
%             i_neq_left=i_neq_left+1;
%         end
%     end
% end
% Bneq(i_neq_right:(i_neq_left-1),1)=0;
% i_neq_right=i_neq_left;

%(17)(19)
for i=1:1:n_path
    for f1=1:1:f_total
        for f2=1:1:f_total
            for a=1:1:n_category
                if ismember(i,A{a})
                    Aneq(i_neq_left,g3(f1,f2,a))=1;
                    Aneq(i_neq_left,g2(f1,f2,i))=-1;
                    Bneq(i_neq_left,1)=0;
                    i_neq_left=i_neq_left+1;                    
                else
                    Aneq(i_neq_left,g3(f1,f2,a))=1;
                    Aneq(i_neq_left,g2(f1,f2,i))=1;
                    Bneq(i_neq_left,1)=1;
                    i_neq_left=i_neq_left+1;                    
                end
            end
        end
    end
end

for a=1:1:n_category
    for f1=1:1:f_total
        for f2=1:1:f_total
            Bneq(i_neq_left,1)=n_path-1;
            Aneq(i_neq_left,g3(f1,f2,a))=-1;
            for i=1:n_path
                if ismember(i,A{a})
                    Aneq(i_neq_left,g2(f1,f2,i))=1;
                else
                    Aneq(i_neq_left,g2(f1,f2,i))=-1;
                    Bneq(i_neq_left,1)=Bneq(i_neq_left,1)-1;
                end
            end
            i_neq_left=i_neq_left+1;
        end
    end
end

i_neq_right=i_neq_left;
for a=1:1:n_category
    for f1=1:1:f_total
        for f2=1:1:f_total
            Aneq(i_neq_left,g3(f1,f2,a))=-1;
        end
    end
    i_neq_left=i_neq_left+1;
end
Bneq(i_neq_right:(i_neq_left-1),1)=-1;

% (9)
for i=1:1:n_path
    Aeq(i_eq_left,x(i,l_max,dest_aug(i)))=1;
    i_eq_left=i_eq_left+1;
end
for i=1:n_source
    for j=1:n_dest
        Aeq(i_eq_left,x((i-1)*n_dest+j,1,i))=1;
        i_eq_left=i_eq_left+1;
    end
end
Beq(i_eq_right:(i_eq_left-1),1)=1;

if length(Aneq(1))<total_v
Aneq(1,total_v)=0;
end


%% minimize order representation
% 
f=sparse(total_v,1);
f(delta_k(1):delta_k(n_dummy),1)=1;
lb=zeros(total_v,1);
ub=ones(total_v,1);

ctype=[];
for i=1:1:total_v
    ctype=strcat(ctype,'B');
end
options=cplexoptimset('MaxIter',1000000000,'MaxNodes',100000000,'MaxTime',10800);
[result,~,exit]=cplexmilp(f,Aneq,Bneq,Aeq,Beq,[],[],[],lb,ub,ctype,[],options);

% tic
% options=optimoptions('intlinprog','MaxTime',36000);
% [result,fval,exitflag,output]=intlinprog(f,1:total_v,Aneq,Bneq,Aeq,Beq,lb,ub,options);
% toc
end
%% minimize size representation
% f=sparse(total_v,1);
% for f1=1:1:f_total
%     for f2=1:1:f_total
%         f(delta_ij(f1,f2):delta_ij(f1,f2),1)=1;
%     end
% end
% lb=zeros(total_v,1);
% ub=ones(total_v,1);
% ctype=[];
% for i=1:1:total_v
%     ctype=strcat(ctype,'B');
% end
% options=cplexoptimset('MaxIter',1000000000,'MaxNodes',100000000,'MaxTime',21600);
% [result,flag,exit]=cplexmilp(f,Aneq,Bneq,Aeq,Beq,[],[],[],lb,ub,ctype,[],options);
% end

function index=x(i,j,f)
    global l_max
    global f_total
    index=(i-1)*l_max*f_total+(j-1)*f_total+f;
end

function index=delta_k(k)
    global n_path
    global l_max
    global f_total
    index=n_path*l_max*f_total+k;
end

function index=delta_ij(f1,f2)
    global n_path
    global l_max
    global f_total
    global n_dummy
    index=n_path*l_max*f_total+n_dummy+(f1-1)*f_total+f2;
end

function index=g1(f1,f2,j1,j2,i)
    global n_path
    global l_max
    global f_total
    global n_dummy
    index=n_path*l_max*f_total+n_dummy+f_total*f_total+ ...,
        (f1-1)*f_total*l_max*l_max*n_path+(f2-1)*l_max*l_max*n_path+ ...,
        (j1-1)*l_max*n_path+(j2-1)*n_path+i;
end

function index=g2(f1,f2,i)
    global n_path
    global l_max
    global f_total
    global n_dummy
    index=n_path*l_max*f_total+n_dummy+f_total*f_total+ ...,
        f_total*f_total*l_max*l_max*n_path+ ...,
        (f1-1)*f_total*n_path+(f2-1)*n_path+i;
end

function index=g3(f1,f2,A)
    global n_path
    global l_max
    global f_total
    global n_dummy
    global n_category
    index=n_path*l_max*f_total+n_dummy+f_total*f_total+ ...,
        f_total*f_total*l_max*l_max*n_path+f_total*f_total*n_path+ ...,
        (f1-1)*f_total*n_category+(f2-1)*n_category+A;
end
    
    