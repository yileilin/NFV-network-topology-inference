function [G_aug_SAP]=sum_SAP(result,weight)
global n_path;
global l_max;
global f_total;

G_aug_SAP=sparse(f_total,f_total);
p = result(1:n_path*l_max*f_total);
p = reshape(p,[f_total l_max n_path]);
paths = zeros(n_path,l_max); % paths(i,j): j-th letter in string i (0 if none)
for i=1:n_path
    for j=1:l_max
        letter = find(p(:,j,i));
        if length(letter) > 1
            disp(['error: i=' num2str(i) ', j=' num2str(j)])
        end
        if ~isempty(letter)
            paths(i,j) = letter;
        end
    end
end
% ignore the non-existing letters to get final augmented strings:
Paths = cell(1,n_path);
for i=1:n_path
    Paths{i} = paths(i,(paths(i,:)>0));
end
size=length(G_aug_SAP);
count=cell(f_total,f_total);
for i=1:n_path
    path=Paths{1,i};
    m=length(path);
    for j=1:m-1
        index_1=path(j);
        index_2=path(j+1);
        count{index_1,index_2}=[count{index_1,index_2},i];
    end
end
for i=1:size
    count{i,i}=[];
end

frequency=zeros(1,2^n_path-1);
matrix=construct_A(n_path);
for i=1:f_total
    for j=1:f_total
        for k=1:2^n_path-1
            if isequal(count{i,j},matrix{1,k})
                frequency(k)=frequency(k)+1;
            end
        end
    end
end

for i=1:f_total
    for j=1:f_total
        for k=1:2^n_path-1
            if isequal(count{i,j},matrix{1,k})
                G_aug_SAP(i,j)=weight(k)/frequency(k);
                if weight(k)==0
                    G_aug_SAP(i,j)=0.00001;
                end
            end
        end
    end
end

for i=f_total:-1:1
    judge_1=find(G_aug_SAP(i,:)>0);
    judge_2=find(G_aug_SAP(:,i)>0);
    if isempty(judge_1) && isempty(judge_2)
        G_aug_SAP(i,:)=[];
        G_aug_SAP(:,i)=[];
    end
end

% SAP_source=1;
% SAP_dest=2:(n_path+1);
% N_VNF=5;
% SAP_VNF=cell(1,N_VNF);
% for i=1:N_VNF
%     SAP_VNF{1,i}=n_path+1+i;
% end
% figure(5)
% G1=digraph(SAP_graph);
% f=plot(G1,'LineWidth',1.5);
% labelnode(f,[1:13],{'' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '' ''});
% labelnode(f,[1:13],{'s' 't1' 't2' 't3' 't4' 't5' 'f1' 'f2' 'f3' 'f4' 'f5' '' ''});
% highlight(f,7:11,'NodeColor','g');
% highlight(f,1,'NodeColor','k');
% highlight(f,2:6,'NodeColor','r');


end

            