function [ H ] = graph_configuration( G)
%if edge weight==0, then edge weight=-1;
%VNF_table records the node which owns that type of VNF
%one VNF does not allow to have replica at the same node

WEIGHT_MAX=5;%weight is the value between 0 and WEIGHT_RANGE
WEIGHT_MIN=1;
p_VNF=0.025;%fraction of nodes that can host VNF
capacity=1; % maximum number of VNF instances per node
n=size(G,1);
H=sparse(n,n);

for i=1:1:n
    for j=1:1:n
        if G(i,j)>0
            H(i,j)=WEIGHT_MIN+(WEIGHT_MAX-WEIGHT_MIN)*rand;
        end
    end
end

% n_list=ceil(p_VNF*n); % #servers
% VNF_list=zeros(capacity+1,n_list);
% %line 1: node indices of servers    line 2-capacity+1: node indices of VNF
% a=randperm(n);
% VNF_list(1,:)=a(1:n_list); % node indices of servers
% i=2;
% j=1;
% for index=1:1:N_VNF
%     VNF_list(i,j)=index;
%     i=i+1;
%     if(i>capacity+1)
%         j=j+1;
%         i=2;
%     end
% end
% 
% for i=(ceil(N_VNF/capacity)+1):1:n_list
%     k=randi(capacity);
%     a=randperm(N_VNF);
%     for j=1:1:k 
%         VNF_list(j+1,i)=a(j);
%     end
% end
% 
% index=ones(1,N_VNF);
% VNF_table=zeros(N_VNF,n_list);
% for i=2:1:capacity+1
%     for j=1:1:n_list
%         if VNF_list(i,j)>0
%             VNF_table(VNF_list(i,j),index(VNF_list(i,j)))=VNF_list(1,j);
%             index(VNF_list(i,j))=index(VNF_list(i,j))+1;
%         end
%     end
% end         

end

