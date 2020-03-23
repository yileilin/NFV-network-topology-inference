AS = '6461';
load(['AS' AS '.mat']);
%% setting
n_path=4;
l_max=11;
n_dummy=5;
n_dest=2;
n_source=2;

n_path=6;
l_max=12;
n_dummy=5;
n_dest=3;
n_source=2;

N_VNF=5;
n_path=3;
n_tree=2;

%% main
G=A;
runs=20;
for int=1:runs
    VNF_set=[];
    dest_set=[];
    n_node=length(G(1,:));
    for i=1:n_node
        if length(find(G(i,:)>0))>6
            VNF_set=[VNF_set,i];
        else
            if length(find(G(i,:)>0))<3
                dest_set=[dest_set,i];
            end
        end
    end

    G = graph_configuration(G);
    VNF_table=zeros(5,5);
    random_array=randperm(length(VNF_set));
    for i=1:1:N_VNF
        VNF_table(i,1)=VNF_set(random_array(i));
    end

    path_table=cell(n_tree,n_path);
    random_array=randperm(length(dest_set));
    source=zeros(n_tree,1);
    dest=zeros(n_tree,n_path);
    for index=1:n_tree
        for i=1:1:n_path
            chain_length=randi(N_VNF);  %chain length from 1 to N_VNF
            sort=randperm(N_VNF);
            path_table{index,i}=sort(1:chain_length);
        end

        source(index,1)=dest_set(random_array(index));
        dest(index,:)=dest_set(random_array((n_tree+1):(n_path+n_tree))); %initial source and distinations
    end

    path_seq=cell(n_tree,n_path);
    for i=1:n_tree
        for i_path=1:1:n_path %find shortest path length and path sequence
            VNF_chain=path_table{i_path};
            n_chain=length(VNF_chain);
            s=source(i,1);
            path_seq{i,i_path}=[];
            path_seq{i,i_path}=[path_seq{i,i_path},s];
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
                path_seq{i,i_path}=[path_seq{i,i_path},sp{shortest_node}];
                s=shortest_node;
            end
            [spcost, sp] = Dijkstra_source(G, s);
            sp{dest(i,i_path)}(1)=[];
            path_seq{i,i_path}=[path_seq{i,i_path},sp{dest(i,i_path)}];
        end
    end

    [GT_graph,GT_source,GT_dest,GT_VNF,GT_path_seq]=contract_multi(path_table,G,path_seq,VNF_table,source,dest);

    path_length=zeros(1,n_tree*n_path);
    for i=1:n_tree*n_path
        path_length(1,i)=routing(GT_graph,GT_path_seq{1,i});
    end

    sharedpath_length=zeros(1,n_tree*n_path*(n_path*n_tree-1)/2);
    index=1;
    n_node=length(GT_graph(1,:));
    for i=1:1:n_path*n_path-1
        for j=i+1:1:n_path*n_tree  %calculate shared portion between path i and j
            path_i=GT_path_seq{i};
            path_j=GT_path_seq{j};
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
                    sharedpath_length(index)=sharedpath_length(index)+min(graph_i(i_loop,j_loop),graph_j(i_loop,j_loop))*GT_graph(i_loop,j_loop);
                end
            end
            index=index+1;
        end
    end

    F=digraph(GT_graph);
    plot(F);
    Inputs{1,int}.GT_graph=GT_graph;
    Inputs{1,int}.GT_source=GT_source;
    Inputs{1,int}.GT_dest=GT_dest;
    Inputs{1,int}.GT_VNF=GT_VNF;
    Inputs{1,int}.GT_path_seq=GT_path_seq;
    Inputs{1,int}.path_length=path_length;
    Inputs{1,int}.sharedpath_length=sharedpath_length;
    Inputs{1,int}.path_table=path_table;
end

function [infer]=routing(G,path)
    infer=0;
    m=length(path);
    for i=1:m-1
        [spcost, ~] = Dijkstra_source(G, path(1,i));
        spcost(1,path(1,i))=0;
        infer=infer+spcost(1,path(i+1));
    end
end