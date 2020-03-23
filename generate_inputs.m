AS = '6461';
load(['AS' AS '.mat']);
global n_path;
global l_max;
global f_total;
global n_dummy;
global n_category;
runs = 20; % maximum number of Monte Carlo runs
Inputs = cell(1,runs); % Inputs{i,r} is the input for r-th run under setting i

%% generate input for n=3:
setting = 1;
N_VNF=5;
n_path=3;
l_max=10;
n_dummy=3;
f_total=1+n_path+N_VNF+n_dummy; %f_total=1+n_dest+N_VNF+n_dummy

for i=1:1:runs
% tic
[ G, VNF_table, path_table, path_seq, path_length, sharedpath_length, source, dest,groundtruth] = shortest_path_hop( A,n_path,N_VNF );
Inputs{setting,i}.G = G;
Inputs{setting,i}.VNF_table = VNF_table;
Inputs{setting,i}.path_table = path_table;
Inputs{setting,i}.path_seq = path_seq;
Inputs{setting,i}.path_length = path_length;
Inputs{setting,i}.sharedpath_length = sharedpath_length;
Inputs{setting,i}.source = source;
Inputs{setting,i}.dest = dest;
Inputs{setting,i}.groundtruth = groundtruth; 
% toc
end
disp(['done for n_path = ' num2str(n_path)])

%% for n=4:
setting = 1;
N_VNF=5;
n_path=4;
l_max=12;
n_dummy=5;
f_total=1+n_path+N_VNF+n_dummy; %f_total=1+n_dest+N_VNF+n_dummy

for i=1:1:runs
% tic
flag=1;
while(flag)
    [ G, VNF_table, path_table, path_seq, path_length, sharedpath_length, source, dest,groundtruth] = shortest_path_hop( A,n_path,N_VNF );
    usage=zeros(1,N_VNF);
    for j=1:n_path
        for k=1:N_VNF
            if ismember(k,path_table{1,j})
                usage(k)=1;
            end
        end
    end
    ones_matrix=ones(1,N_VNF);
    if usage==ones_matrix
        flag=0;
    end
end
Inputs{setting,i}.G = G;
Inputs{setting,i}.VNF_table = VNF_table;
Inputs{setting,i}.path_table = path_table;
Inputs{setting,i}.path_seq = path_seq;
Inputs{setting,i}.path_length = path_length;
Inputs{setting,i}.sharedpath_length = sharedpath_length;
Inputs{setting,i}.source = source;
Inputs{setting,i}.dest = dest;
Inputs{setting,i}.groundtruth = groundtruth; 
end
disp(['done for n_path = ' num2str(n_path)])

%% for n=5:
setting = 1;
N_VNF=5;
n_path=5;
l_max=13;
n_dummy=6;
f_total=1+n_path+N_VNF+n_dummy; %f_total=1+n_dest+N_VNF+n_dummy

for i=1:1:runs
% tic
[ G, VNF_table, path_table, path_seq, path_length, sharedpath_length, source, dest,groundtruth] = shortest_path_hop( A,n_path,N_VNF );
Inputs{setting,i}.G = G;
Inputs{setting,i}.VNF_table = VNF_table;
Inputs{setting,i}.path_table = path_table;
Inputs{setting,i}.path_seq = path_seq;
Inputs{setting,i}.path_length = path_length;
Inputs{setting,i}.sharedpath_length = sharedpath_length;
Inputs{setting,i}.source = source;
Inputs{setting,i}.dest = dest;
Inputs{setting,i}.groundtruth = groundtruth; 
end
disp(['done for n_path = ' num2str(n_path)])

%% for n=6:
setting = 4;
N_VNF=5;
n_path=6;
l_max=13;
n_dummy=6;
f_total=1+n_path+N_VNF+n_dummy; %f_total=1+n_dest+N_VNF+n_dummy

for i=1:1:runs
% tic
[ G, VNF_table, path_table, path_seq, path_length, sharedpath_length, source, dest,groundtruth] = shortest_path_hop( A,n_path,N_VNF );
Inputs{setting,i}.G = G;
Inputs{setting,i}.VNF_table = VNF_table;
Inputs{setting,i}.path_table = path_table;
Inputs{setting,i}.path_seq = path_seq;
Inputs{setting,i}.path_length = path_length;
Inputs{setting,i}.sharedpath_length = sharedpath_length;
Inputs{setting,i}.source = source;
Inputs{setting,i}.dest = dest;
Inputs{setting,i}.groundtruth = groundtruth; 
end
disp(['done for n_path = ' num2str(n_path)])
