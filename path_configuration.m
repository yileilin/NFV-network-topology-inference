function [ path_table ] = path_configuration( N_VNF,n_path)

path_table=cell(1,n_path);
for i=1:1:n_path
    chain_length=randi(N_VNF);  %chain length from 1 to N_VNF
    sort=randperm(N_VNF);
    path_table{i}=sort(1:chain_length);
end
end

