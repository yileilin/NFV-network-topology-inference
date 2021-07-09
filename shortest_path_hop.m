function [ G, VNF_table, path_table, path_seq, path_length, sharedpath_length, source, dest,groundtruth] = shortest_path_hop( G,n_path,N_VNF )

VNF_set=[];

dest_set=[];
n_node=size(G,1);
for i=1:1:n_node
    if length(find(G(i,:)>0))>8
        VNF_set=[VNF_set,i];
    else
        if length(find(G(i,:)>0))<3
            dest_set=[dest_set,i];
        end
    end
end

[ G, VNF_table ] = graph_configuration( G, N_VNF);
VNF_table=zeros(5,5);
random_array=randperm(length(VNF_set));
for i=1:1:5
    for j=1:1:2
        VNF_table(i,j)=VNF_set(random_array(i+(j-1)*5));
    end
end
[ path_table ] = path_configuration( N_VNF,n_path);
n_path=length(path_table);
random_array=randperm(length(dest_set));

% i=[1:n_node];tem(i)=i;
% for i=1:1:n_node
%     if ismember(n_node+1-i,VNF_table)
%         tem(:,n_node+1-i)=[];
%     end
% end
% random_array=randperm(length(tem));
source=dest_set(random_array(1));
dest=dest_set(random_array(2:n_path+1)); %initial source and distinations

path_seq=cell(1,n_path);
path_length=zeros(1,n_path);

for i_path=1:1:n_path %find shortest path length and path sequence
    VNF_chain=path_table{i_path};
    n_chain=length(VNF_chain);
    s=source;
    path_seq{i_path}=[];
    path_seq{i_path}=[path_seq{i_path},s];
    for j_chain=1:1:n_chain
        [spcost, sp] = Dijkstra_source(G, s);
        d=VNF_chain(j_chain);
        shortest_seg=Inf;
        shortest_node=0;
        index=1;
        while(VNF_table(d,index)>0&&index<=N_VNF)
            number_node=VNF_table(d,index);
            if(spcost(number_node)<shortest_seg)
                shortest_seg=spcost(number_node);
                shortest_node=number_node;        
            end
            index=index+1;
        end
        sp{shortest_node}(1)=[];
        path_seq{i_path}=[path_seq{i_path},sp{shortest_node}];
        path_length(1,i_path)=path_length(1,i_path)+shortest_seg;
        s=shortest_node;
    end
    [spcost, sp] = Dijkstra_source(G, s);
    sp{dest(i_path)}(1)=[];
    path_seq{i_path}=[path_seq{i_path},sp{dest(i_path)}];
    path_length(1,i_path)=path_length(1,i_path)+spcost(dest(i_path));
end

sharedpath_length=zeros(1,n_path*(n_path-1)/2);

index=1;
for i=1:1:n_path-1
    for j=i+1:1:n_path  %calculate shared portion between path i and j
        path_i=path_seq{i};
        path_j=path_seq{j};
        length_i=length(path_i);
        length_j=length(path_j);
        
        graph_i=sparse(n_node,n_node);
        graph_j=sparse(n_node,n_node);
        for i_loop=1:1:length(path_i)-1
            graph_i(path_i(i_loop),path_i(i_loop+1))=graph_i(path_i(i_loop),path_i(i_loop+1))+1;
        end
        for j_loop=1:1:(length_j-1)
            graph_j(path_j(j_loop),path_j(j_loop+1))=graph_j(path_j(j_loop),path_j(j_loop+1))+1;
        end
        for i_loop=1:1:n_node
            for j_loop=1:1:n_node
                sharedpath_length(index)=sharedpath_length(index)+min(graph_i(i_loop,j_loop),graph_j(i_loop,j_loop))*G(i_loop,j_loop);
            end
        end
        index=index+1;
    end
end

groundtruth=sparse(n_node,n_node); %condense the graph into groundtruth
% note: should contract edge sequences into edges
for i_path=1:1:n_path
    path_i=path_seq{i_path};
    length_i=length(path_i);
    for index=1:1:length_i-1
        groundtruth(path_i(index),path_i(index+1))=G(path_i(index),path_i(index+1));
    end
end

% groundtruth3=groundtruth;
% i=1;
% index=1;
% while index<=length(groundtruth)
%     b=find(groundtruth3(i,:)>1e-8);
%     c=find(groundtruth3(:,i)>1e-8);
%     if isempty(b) && isempty(c)
%         groundtruth3(i,:)=[];
%         groundtruth3(:,i)=[];
%         i=i-1;
%     end
%     i=i+1;
%     index=index+1;
% end
% 
% 
% G1=digraph(groundtruth3);
% plot(G1);
% groundtruth2=groundtruth3;
% for i=1:1:length(groundtruth3)
%     temp1=find(groundtruth3(:,i)>1e-8);
%     if(length(temp1)==1)
%         temp2=find(groundtruth3(i,:)>1e-8);
%         if(length(temp2)==1&&temp1(1)~=temp2(1))
%             groundtruth3(temp1,temp2)=groundtruth3(temp1(1),i)+groundtruth3(i,temp2(1));
%             groundtruth3(temp1(1),i)=0;
%             groundtruth3(i,temp2(1))=0;
%         end
%     end
% end
% 
% G1=digraph(groundtruth3);
% plot(G1);
% 
% %groundtruth3=groundtruth2;
% i=1;
% index=1;
% while index<=length(groundtruth3);
%     b=find(groundtruth3(i,:)>1e-8);
%     c=find(groundtruth3(:,i)>1e-8);
%     if isempty(b) && isempty(c)
%         groundtruth3(i,:)=[];
%         groundtruth3(:,i)=[];
%         i=i-1;
%     end
%     i=i+1;
%     index=index+1;
% end
% 
% G1=digraph(groundtruth3);
% plot(G1);
% 
% need=zeros(1,n_node);
% for i=1:1:n_node
%     for j=1:1:n_node
%         if groundtruth2(i,j)>1e-8
%             need(1,i)=1;
%             need(1,j)=1;
%         end
%     end
% end
% 
% hop=0;
% for i=1:1:n_path
%     p=path_seq{i};
%     for j=1:1:length(p)
%         if need(p(j)) 
%             hop=hop+1;
%         end
%     end
% end
% hop=hop-n_path;
% 
% n_vertice=length(find(need==1));
% 
% n_edge=nnz(groundtruth2);
% 
% total_weight=sum(groundtruth(find(groundtruth>1e-8)),1);
end

