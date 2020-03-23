clear all
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64')

global n_path;
global l_max;
global f_total;
global n_dummy;
global n_category;

setting = 1;
N_VNF=5;
n_path=6;
l_max=18;
n_dummy=5;
f_total=1+n_path+N_VNF+n_dummy; %f_total=1+n_dest+N_VNF+n_dummy

load('./multi_round1/Inputs.mat'); % Inputs
load('common_6.mat');

Outputs=cell(5,20);
sequence=cell(1,20);

blank=cell(1,n_path);

% a structure, 1:groundtruth  2:RNJ++  3:CE  4:CE++  5.SAP
for i=1:20
    G = Inputs{setting,i}.G;
    VNF_table = Inputs{setting,i}.VNF_table;
    path_table = Inputs{setting,i}.path_table;
    path_seq = Inputs{setting,i}.path_seq;
    path_length = Inputs{setting,i}.path_length;
    sharedpath_length = Inputs{setting,i}.sharedpath_length;
    source = Inputs{setting,i}.source;
    dest = Inputs{setting,i}.dest;
    groundtruth = Inputs{setting,i}.groundtruth;
    
    [weight,fval,exitflag] = min_weight_assignment_LP(path_length,sharedpath_length);
    %weight = category_classify( G , path_seq );
    %GT
    [GT_graph,GT_source,GT_dest,GT_VNF]=contract(path_table,groundtruth,path_seq,VNF_table,source,dest,1);
    Outputs{1,i}.adjacency=GT_graph;
    Outputs{1,i}.source=GT_source;
    Outputs{1,i}.destination=GT_dest;
    Outputs{1,i}.vnf=GT_VNF;
    
    %RNJ++
    [RNJ_tree,RNJ_source,RNJ_dest]=RNJ_plus(path_length,sharedpath_length,dest,path_table);
    Outputs{2,i}.adjacency=RNJ_tree;
    Outputs{2,i}.source=RNJ_source;
    Outputs{2,i}.destination=RNJ_dest;
    Outputs{2,i}.vnf=blank;
    
    %CE
    [CE_graph,CE_source,CE_dest]=sum_CE(weight,path_table);
    Outputs{3,i}.adjacency=CE_graph;
    Outputs{3,i}.source=CE_source;
    Outputs{3,i}.destination=CE_dest;
    Outputs{3,i}.vnf=blank;

    %CE++
    [CE_plus_graph,CE_plus_source,CE_plus_dest,CE_plus_VNF]=sum_CE_plus(weight,path_table);
    Outputs{4,i}.adjacency=CE_plus_graph;
    Outputs{4,i}.source=CE_plus_source;
    Outputs{4,i}.destination=CE_plus_dest;
    Outputs{4,i}.vnf=CE_plus_VNF;
    
    %SAP
    sequence{1,i}=[];
    [Aeq,Beq,Aneq,Bneq,result,exit]=SAP(weight,path_table,n_path,N_VNF,dest,n_dummy,l_max,Aneq1,Bneq1,i_neq_left);
    sequence{1,i}=result;
    [SAP_graph,SAP_source,SAP_dest,SAP_VNF]=sum_SAP(result);
 %   [SAP_graph,SAP_source,SAP_dest,SAP_VNF]=sum_SAP(sequence{1,i});
    Outputs{5,i}.adjacency=SAP_graph;
    Outputs{5,i}.source=SAP_source;
    Outputs{5,i}.destination=SAP_dest;
    Outputs{5,i}.vnf=SAP_VNF;
    %[ val(1,i) ] = graph_edit_distance( Outputs{1,i}, Outputs{5,i} )
    save('./round2/order10800.mat','Outputs','sequence');
end
