load('./multi_round1/order_2_3.mat')
load('./multi_round1/Inputs_2_3.mat')
sap=[];
gt=[];
ce=[];
rnj=[];

for index=1:5
    %GT
    n=length(Outputs{1,index}.adjacency);
    array_gt=zeros(1,n);
    k=1;
    for i=1:n
        for j=1:n
            if (Outputs{1,index}.adjacency(i,j)>0)
                array_gt(1,k)=Outputs{1,index}.adjacency(i,j);
                k=k+1;
            end
        end
    end
    [gt]=[gt,array_gt];
    
    %RNJ
    n=length(Outputs{2,index}.adjacency);
    array_rnj=zeros(1,n);
    k=1;
    for i=1:n
        for j=1:n
            if (Outputs{2,index}.adjacency(i,j)>0)
                array_rnj(1,k)=Outputs{2,index}.adjacency(i,j);
                k=k+1;
            end
        end
    end
    [rnj]=[rnj,array_rnj];
    
    %CE
    n=length(Outputs{3,index}.adjacency);
    array_ce=zeros(1,n);
    k=1;
    for i=1:n
        for j=1:n
            if (Outputs{3,index}.adjacency(i,j)>0)
                array_ce(1,k)=Outputs{3,index}.adjacency(i,j);
                k=k+1;
            end
            if (Outputs{3,index}.adjacency(i,j)<0)
                array_ce(1,k)=0.000001;
                k=k+1;
            end
        end
    end    
    [ce]=[ce,array_ce];
    
    %SAP
    n=length(Outputs{4,index}.adjacency(1,:));
    array_sap=zeros(1,n);
    k=1;
    for i=1:n
        for j=1:n
            if (Outputs{4,index}.adjacency(i,j)>0)
                array_sap(1,k)=Outputs{4,index}.adjacency(i,j);
                k=k+1;
            end
        end
    end
    [sap]=[sap,array_sap];
    
    
end

figure(1)
A=cdfplot(gt);
hold on
B=cdfplot(rnj);
hold on
C=cdfplot(ce);
hold on
D=cdfplot(sap);
set(A,'LineWidth',1.5);
set(B,'LineWidth',1.5);
set(C,'LineWidth',1.5);
set(D,'LineWidth',1.5);
xlabel('edge weight','FontSize',20);
legend('ground truth','RNJ','CE','SAP','FontSize',12);
set(gca,'LooseInset',get(gca,'TightInset'));

% figure(2)
% %x=['GT','RNJ','CE','SAP'];
% y=[mean(gt),mean(rnj),mean(ce),mean(sap)];
% f1=bar(y);
% 
% figure(3)
% %x=['GT','RNJ','CE','SAP'];
% y=[max(gt),max(rnj),max(ce),max(sap)];
% bar(y);
% 
% figure(4)
% x=[3,4,5];
% y=[2.7637,12.5287,4.746,4.1994;3.3341,13.4767,4.2371,5.2937;2.6168,8.6506,3.2501,3.6073];
% bar(x,y);
% legend('Ground truth','RNJ','CE','SAP');
% xlabel('number of path');
% ylabel('average edge weight');
% 
% figure(5)
% x=[3,4,5];
% y=[3.5333,0.7173,0.5195;3.0421,0.2708,0.5877;2.3058,0.242,0.3785];
% bar(x,y);
% legend('RNJ','CE','SAP')
% xlabel('number of path');
% ylabel('relative error');