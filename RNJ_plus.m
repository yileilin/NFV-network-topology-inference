function [RNJ_tree,RNJ_source,RNJ_dest]=RNJ_plus(path_length,sharedpath_length)
n_path=length(path_length);
% path_length=[14.8871,15.4898,7.0997];
% sharedpath_length=[1.2657,7.0997,1.2657];
delta=0;
RNJ_tree=sparse(n_path,n_path);% add s in the end
flag=1;
D=ones(1,n_path);%1 in, 0 out
index_node=n_path;% indices of f

if n_path==2
    RNJ_tree=zeros(4,4);
    RNJ_source=4;
    RNJ_dest=[1,2];
    RNJ_tree(4,3)=sharedpath_length;
    RNJ_tree(3,1)=path_length(1)-sharedpath_length;
    RNJ_tree(3,2)=path_length(2)-sharedpath_length;
else
    rho=zeros(3,n_path*(n_path-1)/2);
    index=1;
    for i=1:1:n_path-1
        for j=i+1:1:n_path
            rho(1,index)=i;
            rho(2,index)=j;
            rho(3,index)=sharedpath_length(index);
            index=index+1;
        end
    end

    while flag
        max_value=-1;
        position=0;
        for i=1:1:length(rho(1,:))
            if D(rho(1,i))==1 && D(rho(2,i))==1
                if rho(3,i)>max_value
                    max_value=rho(3,i);
                    position=i;
                end
            end    
        end
        index_node=index_node+1;
        i_star=rho(1,position);
        j_star=rho(2,position);
        D(i_star)=0;
        D(j_star)=0;
        RNJ_tree(index_node,i_star)=path_length(i_star)-max_value;
        RNJ_tree(index_node,j_star)=path_length(j_star)-max_value;

        for i=1:1:length(D)
            if rho(1,i)==i_star && D(rho(2,i))==1 && max_value-rho(3,i)<=delta
                D(rho(2,i))=0;
                RNJ_tree(index_node,rho(2,i))=path_length(rho(2,i))-max_value;
            end
        end


        for k=1:1:length(D)
            if D(k)==1
                sum=0;
                for j=1:1:length(rho(1,:))
                    if rho(1,j)==k && rho(2,j)==i_star || rho(1,j)==i_star && rho(2,j)==k || ...,
                        rho(1,j)==k && rho(2,j)==j_star || rho(1,j)==j_star && rho(2,j)==k
                        sum=sum+rho(3,j);
                    end
                end
                rho=[rho,[k;index_node;sum/2]];
            end      
        end
        D(index_node)=1;
        path_length(index_node)=max_value;
        if length(find(D==1))<=1
            RNJ_tree(index_node+1,index_node)=path_length(length(path_length));
            index_node=index_node+1;
            flag=0;
        end
    end
    RNJ_tree(1,index_node)=0;

    RNJ_source=length(RNJ_tree);
    RNJ_dest=[1:n_path];
end

   
end