% load('./data/simulation_test_2.mat');
% setting = 1;
% hop=zeros(1,20);
% n_edge=zeros(1,20);
% n_vertice=zeros(1,20);
% total_weight=zeros(1,20);
% for c=1:1:1
% %[ G, VNF_table, path_table, path_seq, path_length, sharedpath_length, source, dest,groundtruth] = shortest_path_hop( A,n_path,N_VNF );
% %load pre-generated inputs:
% G = Inputs{setting,c}.G;
% VNF_table = Inputs{setting,c}.VNF_table;
% path_table = Inputs{setting,c}.path_table;
% path_seq = Inputs{setting,c}.path_seq;
% path_length = Inputs{setting,c}.path_length;
% sharedpath_length = Inputs{setting,c}.sharedpath_length;
% source = Inputs{setting,c}.source;
% dest = Inputs{setting,c}.dest;
% groundtruth = Inputs{setting,c}.groundtruth;
function [groundtruth3,GT_source,GT_dest,GT_VNF]=contract(path_table,groundtruth,path_seq,VNF_table,source,dest,dup)
n_path=length(path_table);
N_VNF=length(VNF_table(:,1));
%server
n_node=length(groundtruth);
co_node=zeros(1,n_node);
for i=1:n_node
    co_node(i)=i;
end

index=1;
server=[];
for h=1:1:n_path
    m=length(path_seq{h});
    n=length(path_table{h});
    p=path_seq{h};
    j=1;
    k=1;
    while j<=n
        while ~ismember(p(k),VNF_table(path_table{h}(j),:))&&~ismember(p(k),dest)&&~ismember(p(k),source)
            k=k+1;
        end
        server(index)=p(k);
        index=index+1;
        j=j+1;
    end
end

for i=1:1:length(groundtruth)
    temp1=find(groundtruth(:,i)>1e-8);
    if(length(temp1)==1)
        temp2=find(groundtruth(i,:)>1e-8);
        if(length(temp2)==1&&temp1(1)~=temp2(1)&&~ismember(i,server)&&~ismember(p(k),dest)&&~ismember(p(k),source))%server
            groundtruth(temp1,temp2)=groundtruth(temp1(1),i)+groundtruth(i,temp2(1));
            groundtruth(temp1(1),i)=0;
            groundtruth(i,temp2(1))=0;
%             for z=(i+1):n_node
%                 co_node(z)=co_node(z)-1;
%             end
        end
    else
        if(length(temp1)==2)
            temp2=find(groundtruth(i,:)>1e-8);
            if length(temp2)==2 && ~ismember(i,server) && temp1(1)==temp2(1) && temp1(2)==temp2(2)&& ~ismember(i,dest)&& ~ismember(i,source)%server
                groundtruth(temp1(1),temp1(2))=groundtruth(temp1(1),i)+groundtruth(i,temp1(2));
                groundtruth(temp1(2),temp1(1))=groundtruth(temp1(2),i)+groundtruth(i,temp1(1));
                groundtruth(temp1(1),i)=0;
                groundtruth(i,temp1(2))=0;
%                 for z=(i+1):n_node
%                     co_node(z)=co_node(z)-1;
%                 end
            end
        end
    end
            
end

% G1=digraph(groundtruth);
% plot(G1);

groundtruth3=groundtruth;
i=1;
index=1;
while index<=length(groundtruth)
    b=find(groundtruth3(i,:)>1e-8);
    d=find(groundtruth3(:,i)>1e-8);
    if isempty(b) && isempty(d)
        groundtruth3(i,:)=[];
        groundtruth3(:,i)=[];
        for z=(index+1):n_node
            co_node(z)=co_node(z)-1;
        end
        co_node(index)=0;
        i=i-1;
    end
    i=i+1;
    index=index+1;
end

GT_source=co_node(source);
GT_dest=co_node(dest);
GT_VNF=cell(1,N_VNF);
for i=1:N_VNF
    GT_VNF{i}=[];
    for j=1:dup
        GT_VNF{i}=[GT_VNF{i},co_node(VNF_table(i,j))];
    end
end
% N_VNF=5;
% color_VNF=zeros(1,n_path);
% for i=1:1:N_VNF
%     color_VNF(1,i)=co_node(VNF_table(i,1));
% end
% figure(1)
% G1=digraph(groundtruth3);
% f=plot(G1,'LineWidth',1.5);
% %labelnode(f,[1:16],{'f2' '' '' 'f3' 'f4' '' 'f1' 'f5' '' 't1' 't2' 't5' 't4' '' 's' 't3'});
% % labelnode(f,[1:16],{'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''});
% 
% highlight(f,color_VNF,'NodeColor','g');
% highlight(f,co_node(source),'NodeColor','k');
% highlight(f,co_node(dest(1:n_path)),'NodeColor','r');

end