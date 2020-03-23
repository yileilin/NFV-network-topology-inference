function [CE_plus,CE_plus_source,CE_plus_dest,CE_plus_VNF]=sum_CE_plus(weight,path_table)
global n_path;
N_VNF=5;
CE_plus=sparse(N_VNF+n_path+1,N_VNF+n_path+1);
I=find(weight>0.0001);
[path_index]=construct_A(n_path);
flag=zeros(1,2^n_path-1);
flag(I)=1;

edge_cate=cell(N_VNF+2+n_path,N_VNF+2+n_path);

%CE_plus   1:source  2-path+1:dest   VNF=VNF+1+n_path
%dummy:from 1+n_path+N_VNF+1
CE_plus(1,2+n_path+N_VNF)=1;
for j=1:1:n_path
    p=path_table{1,j};
    length_p=length(p);
    CE_plus(2+n_path+N_VNF,p(1)+1+n_path)=1;
    edge_cate{2+n_path+N_VNF,p(1)+1+n_path}=[edge_cate{2+n_path+N_VNF,p(1)+1+n_path},j];
    if length_p>1
        for k=1:1:length_p-1
            CE_plus(p(k)+1+n_path,p(k+1)+1+n_path)=1;
            edge_cate{p(k)+1+n_path,p(k+1)+1+n_path}=[edge_cate{p(k)+1+n_path,p(k+1)+1+n_path},j];
        end
    end
    CE_plus(p(length_p)+1+n_path,j+1)=1;
    edge_cate{p(length_p)+1+n_path,j+1}=[edge_cate{p(length_p)+1+n_path,j+1},j];
end
flag(1:n_path)=0;
flag(2^n_path-1)=0;
for j=1:1:2+n_path+N_VNF
    for k=1:1:2+n_path+N_VNF
        if edge_cate{j,k}
            for l=1:1:2^n_path-1
                if isequal(path_index{1,l},edge_cate{j,k})
                    flag(l)=0;
                end
            end
        end
    end
end
n_category=length(find(flag>0));

    
    
if(n_category>0)
    clique_node=ceil(sqrt(n_category));
if(clique_node*(clique_node-1)<n_category)
    clique_node=clique_node+1;
end

CE_graph=sparse(clique_node+1,clique_node+1);
random_array=randperm(n_category);

record_sequence=zeros(1,n_category);%i-th records what category is going to assign
index=1;
for i=1:1:2^n_path-1
    if flag(i)==1
        record_sequence(1,random_array(index))=i;
        index=index+1;
    end
end

CE_category=cell(clique_node+1,clique_node+1);
i=1;
j=2;
for k=1:1:n_category
    CE_graph(i,j)=weight(record_sequence(k));
    CE_category{i,j}=path_index{record_sequence(k)};
    j=j+1;
    if(j==i)
        j=j+1;
    end
    if(j==clique_node+1)
        i=i+1;
        j=1;
    end
end

for k=1:1:n_path
    %extract all the edges in the clique that traverse path k
    temp=sparse(clique_node,clique_node);
    for i=1:1:clique_node
        for j=1:1:clique_node
            if length(CE_category{i,j})>0
                if ismember(k,CE_category{i,j})
                    temp(i,j)=1;
                end
            end
        end
    end
    flag=1;
    
    while(flag)
        I=find(temp==1);
        if ~isempty(I)
        for i=1:1:length(I)
            flag2=1;
            for j=1:1:length(I)
                if(i~=j)
                    n1=I(i);
                    n2=I(j);
                    p1=floor((n1-1)/clique_node)+1;
                    q1=mod(n1,clique_node);
                    if(q1==0)
                        q1=clique_node;
                    end
                    p2=floor((n2-1)/clique_node)+1;
                    q2=mod(n2,clique_node);
                    if(q2==0)
                        q2=clique_node;
                    end
                    if(p1==q2)
                        temp(q1,p1)=0;
                        temp(q2,p2)=0;
                        temp(q1,p2)=1;
                        flag2=0;
                        break;
                    end
                end
                if(i==length(I)&&j==length(I))
                    flag=0;
                end
            end
            if(flag2==0)
                break;
            end
                
        end
        else
            flag=0;
        end
    end
    I=find(temp==1);
    for i=1:1:length(I)
        p=floor((I(i)-1)/clique_node)+1;
        q=mod(I(i),clique_node);
        if(q==0)
            q=clique_node;
        end
        CE_graph(clique_node+1,q)=1;%minus
        CE_graph(p,clique_node+1)=1;%minus
        CE_category{clique_node+1,q}=[CE_category{clique_node+1,q},k];
        CE_category{p,clique_node+1}=[CE_category{p,clique_node+1},k];
    end
end

%N_VNF+2+n_path
for i=1:1:clique_node
    for j=1:1:clique_node
        CE_plus(N_VNF+2+n_path+i,N_VNF+2+n_path+j)=CE_graph(i,j);
    end
end

for i=1:1:clique_node
    CE_plus(N_VNF+2+n_path,i+N_VNF+2+n_path)=CE_graph(clique_node+1,i);
    CE_plus(i+N_VNF+2+n_path,N_VNF+2+n_path)=CE_graph(i,clique_node+1);
end

CE_plus_source=1;
CE_plus_dest=2:(n_path+1);
CE_plus_VNF=cell(1,n_path);
for i=1:N_VNF
   CE_plus_VNF{1,i}=1+n_path+i; 
end


% figure(3)
% G1=digraph(CE_plus);
% f=plot(G1,'LineWidth',1.5);
% 
% labelnode(f,[1:16],{' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '' '' '' '' ''});
% labelnode(f,[1:16],{'s' 't1' 't2' 't3' 't4' 't5' 'f1' 'f2' 'f3' 'f4' 'f5' '' '' '' '' ''});
% 
% highlight(f,(2+n_path:1+n_path+N_VNF),'NodeColor','g');
% highlight(f,1,'NodeColor','k');
% highlight(f,2:1+n_path,'NodeColor','r');
end