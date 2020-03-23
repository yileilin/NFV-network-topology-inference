%configuration
function [CE_graph,CE_source,CE_dest]=sum_CE(weight)
global n_source
global n_dest
global n_path
I=find(weight>0.0001);
n_category=length(I);
clique_node=ceil(sqrt(n_category));
if(clique_node*(clique_node-1)<n_category)
    clique_node=clique_node+1;
end


[path_index]=construct_A(n_path);
flag=zeros(1,2^n_path-1);
flag(I)=1;
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
    end
    I=find(temp==1);
    for i=1:1:length(I)
        p=floor((I(i)-1)/clique_node)+1;
        q=mod(I(i),clique_node);
        if(q==0)
            q=clique_node;
        end
        CE_graph(clique_node+1,q)=0.00001;
        CE_graph(p,clique_node+1)=0.00001;
        CE_category{clique_node+1,q}=[CE_category{clique_node+1,q},k];
        CE_category{p,clique_node+1}=[CE_category{p,clique_node+1},k];
    end
end
for i=1:n_source
    CE_graph(clique_node+1+i,clique_node+1)=0.00001;
end
for i=1:1:n_dest
    CE_graph(clique_node+1,clique_node+1+n_source+i)=0.00001;
end
CE_graph(clique_node+1+n_source+n_dest,1)=0;
CE_source=(clique_node+2):(clique_node+n_source+1);
CE_dest=(clique_node+n_source+2):(clique_node+1+n_source+n_dest);



% figure(2)
% G1=digraph(CE_graph);
% f=plot(G1,'LineWidth',1.5);
% labelnode(f,[1:clique_node+2+n_path],{'' '' '' '' '' '' 's' 't1' 't2' 't3' 't4' 't5'})
% %labelnode(f,[1:clique_node+2+n_path],{'' '' '' '' '' '' '' '' '' '' '' ''});
% 
% highlight(f,clique_node+2,'NodeColor','k');
% highlight(f,(clique_node+3):(clique_node+2+n_path),'NodeColor','r');
end

