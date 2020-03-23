function [groundtruth,GT_source,GT_dest,GT_VNF,GT_path_seq]=contract_multi_new(path_table,G,path_seq,VNF_table,source,dest)
groundtruth=sparse(length(G(1,:)),length(G(1,:)));
n_tree=length(dest(:,1));
n_dest=length(dest(1,:));
GT_path_seq=cell(1,n_tree*n_dest);
co_node=zeros(1,length(G(1,:)));
for i=1:length(source)
    co_node(1,source(i,1))=1;
end
for i=1:length(dest(1,:))
    co_node(1,dest(1,i))=1;
end
for i=1:length(VNF_table(:,1))
    co_node(1,VNF_table(i,1))=1;
end

for i=1:length(G(1,:))
    if co_node(1,i)==0
        node_list=[];
        for i_path=1:n_tree
            for j_path=1:n_dest
                path=path_seq{i_path,j_path};
                if ismember(i,path)
                    X=find(path==i);
                    for j=1:length(X)
                        node_list=[node_list,path(X(j)-1)];
                        node_list=[node_list,path(X(j)+1)];
                    end   
                end
            end
        end
        if length(unique(node_list))>2
            co_node(1,i)=1;
        end
    end
end

for i_path=1:n_tree
    for j_path=1:n_dest
        path=path_seq{i_path,j_path};
        i=1;
        while(i<length(path))
            sum=0;
            j=i+1;
            sum=sum+G(path(i),path(j));
            while co_node(path(j))==0
                sum=sum+G(path(j),path(j+1));
                j=j+1;
            end
            groundtruth(path(i),path(j))=sum;
            i=j;
        end
    end
end

n_node=length(G(1,:));
co_node=zeros(1,n_node);
for i=1:n_node
    co_node(i)=i;
end
i=1;
index=1;
while index<=length(G)
    b=find(groundtruth(i,:)>1e-8);
    d=find(groundtruth(:,i)>1e-8);
    if isempty(b) && isempty(d)
        groundtruth(i,:)=[];
        groundtruth(:,i)=[];
        for z=(index+1):n_node
            co_node(z)=co_node(z)-1;
        end
        co_node(index)=0;
        i=i-1;
    end
    i=i+1;
    index=index+1;
end
N_VNF=5;
GT_source=co_node(source(:,1));
GT_dest=co_node(dest(1,:));
GT_VNF=cell(1,N_VNF);
for i=1:N_VNF
    GT_VNF{i}=[];
    for j=1:1
        GT_VNF{i}=[GT_VNF{i},co_node(VNF_table(i,j))];
    end
end
for i=1:n_tree
    for j=1:n_dest
        path=path_seq{i,j};
        GT_path_seq{1,(i-1)*n_dest+j}=[];
        for i_length=1:length(path)
            if co_node(path(i_length))~=0
                GT_path_seq{1,(i-1)*n_dest+j}=[GT_path_seq{1,(i-1)*n_dest+j},co_node(path(i_length))];
            end
        end
    end
end