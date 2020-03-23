
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64')
global n_path;
global l_max;
global f_total;
global n_dummy;
global n_dest;
global n_category;
global n_source

N_VNF=5;
n_path=2;
l_max=9;
n_dummy=2;
n_dest=1;
n_source=2;
f_total=n_source+n_dest+N_VNF+n_dummy; %f_total=1+n_dest+N_VNF+n_dummy

load('./multi_round_20_runs/Inputs_2_1.mat'); % Inputs
load('common_2_1.mat');
%load('./multi_round1/order_2_3.mat');
%dup=1;
% Outputs=cell(4,20);
% sequence=cell(1,20);

blank=cell(1,n_path);
% a structure, 1:groundtruth  2:RNJ++  3:CE  4:CE++  5.SAP
for i=1:20
    GT_graph = Inputs{1,i}.GT_graph;
    GT_source = Inputs{1,i}.GT_source;
    GT_dest=Inputs{1,i}.GT_dest;
    GT_VNF=Inputs{1,i}.GT_VNF;
    path_seq=Inputs{1,i}.GT_path_seq;
    path_length = Inputs{1,i}.path_length;
    sharedpath_length = Inputs{1,i}.sharedpath_length;
    path_table=Inputs{1,i}.path_table;

% H=digraph(groundtruth);
% plot(H);
    [weight,~,~] = min_weight_assignment_LP(path_length,sharedpath_length);
%     %weight = category_classify( GT_graph , path_seq );
    Outputs{1,i}.adjacency=GT_graph;
    Outputs{1,i}.source=GT_source;
    Outputs{1,i}.destination=GT_dest;
    Outputs{1,i}.vnf=GT_VNF;
%     
    %REA
%     [tree,tree_weight,R,pathlength,sharedlength,s2b2]=REA_pre(path_length,sharedpath_length);
%     [ G, G_weight ] = REA_plus( tree, tree_weight, R, pathlength, sharedlength, s2b2 );
%     Outputs{2,i}.adjacency=G_weight;
%     Outputs{2,i}.source=[2*n_dest,2*n_dest+1];
%     Outputs{2,i}.destination=[1:n_dest];
%     Outputs{2,i}.vnf=blank;
    if sharedpath_length>0
        G_weight=zeros(4,4);
        G_weight(4,1)=sharedpath_length;
        G_weight(2,4)=path_length(1)-sharedpath_length;
        G_weight(3,4)=path_length(2)-sharedpath_length;
    else
        G_weight=zeros(3,3);
        G_weight(2,1)=path_length(1);
        G_weight(3,1)=path_length(2);
    end
    Outputs{2,i}.adjacency=G_weight;
    Outputs{2,i}.source=[2,3];
    Outputs{2,i}.destination=1;
    Outputs{2,i}.vnf=blank;
%     
%     %CE
    [CE_graph,CE_source,CE_dest]=sum_CE(weight);
    Outputs{3,i}.adjacency=CE_graph;
    Outputs{3,i}.source=CE_source;
    Outputs{3,i}.destination=CE_dest;
    Outputs{3,i}.vnf=blank;
%     
%     %SAP
disp('start');
    sequence{1,i}=[];
     [~,~,~,~,result,exit]=SAP_multi(weight,path_table,Aneq,Bneq,i_neq_left);
     sequence{1,i}=result;
    [SAP_graph]=sum_SAP(sequence{1,i},weight);
 %   [SAP_graph,SAP_source,SAP_dest,SAP_VNF]=sum_SAP(sequence{1,i});
    Outputs{4,i}.adjacency=SAP_graph;
    Outputs{4,i}.source=1:n_source;
    Outputs{4,i}.destination=(n_source+1):(n_source+n_dest);
    Outputs{4,i}.vnf=cell(1,N_VNF);
    for j=1:N_VNF
        Outputs{4,i}.vnf{1,j}=n_source+n_dest+j;
    end
    save('./multi_round_20_runs/order_2_1.mat','Outputs','sequence');
end
