function [ G, G_weight ] = REA_plus( tree, weight, R, pathlength, sharedlength, s2b2 )
% implementing a variation of the REA algorithm (Alg 2) in Sattari et al,
% "Active Learning of Multiple Source Multiple Destination Topologies",
% IEEE Trans. SP, 2014.
% Note: original alg requires a test of quartet type, and only output an
% unweighted graph; we directly uses the measured path lengths and shared
% path lengths, and output a weighted graph that encodes both the topology
% and the link performances. 
% Input:
% tree: adjacency matrix of a directed tree from s1 to all receivers R; tree(i,j)=1 means i->j edge
% exists
% weight: weight(i,j) is the weight of edge (i,j)
% R: node indices of receivers in tree
% pathlength: pathlength(i) is the length of s2-to-R(i) path
% sharedlength: sharedlength(i) is the shared length between path
% s1-to-R(i) and path s2-to-R(i)
% s2b2: length of link (s2,b2), i.e., minimum shared length between s2 and
% any two receivers in R
G = zeros(length(tree) + 2 + length(R)); % original nodes + s2 + b2 + at most one joining point per receiver
G_weight = G;
% initially, G = tree
G(1:length(tree),1:length(tree)) = tree;
G_weight(1:length(tree),1:length(tree)) = weight; 
ispending = zeros(1,length(tree)); % ispending(i) = 1 indicates node i's joining point from s1 and s2 has not been identified
ispending(R) = 1;
actual = zeros(1,length(tree)); % receiver index (in 1,...,length(R)) each pending node represents
actual(R) = 1:length(R); 
dist = zeros(1,length(tree)); % distance from a pending node to the receiver it represents
% assume in G, node 1,...,length(tree) are nodes in tree, node
% length(tree)+1 is s2, node length(tree)+2 is b2.
s2 = length(tree)+1; 
b2 = length(tree)+2;
G(s2,b2) = 1;
if s2b2 < min(pathlength - sharedlength) %b2 is above the first joining point    
    G_weight(s2,b2) = s2b2;
end % else: G_weight(s2,b2) = 0 (b2 will be contracted with s2 later)
j_new = length(tree) + 2 + 1; % node index of the next joining point in G
loop=0;
while any(ispending)
    for i=1:length(tree)
        if any(tree(i,:)) && all(tree(i,:) - ispending <= 0) % i has children and all children of node i are pending nodes
            P = i; break;
        end
    end
    i_above = 0; % index of the receiver represented by one of P's children, for which the joining point is above P
%     d_above: distance between P and that receiver
    children = find(tree(P,:)); 
    for Ri = children
        if sharedlength(actual(Ri)) >= dist(Ri) + weight(P,Ri) % Ji is above P
            i_above = actual(Ri); d_above = dist(Ri) + weight(P,Ri);
        else % Ji is below P, i.e., on link (P, Ri)
            G(P,j_new) = 1; 
            G(j_new,Ri) = 1;            
            G(b2,j_new) = 1;
            G(P,Ri) = 0; 
            G_weight(j_new,Ri) = sharedlength(actual(Ri)) - dist(Ri);
            G_weight(P,j_new) = G_weight(P,Ri) - G_weight(j_new,Ri); 
            G_weight(b2,j_new) = pathlength(actual(Ri)) - G_weight(s2,b2) - sharedlength(actual(Ri)); 
            G_weight(P,Ri) = 0;
            j_new = j_new + 1;
        end
    end
    ispending(children) = 0; 
    if i_above > 0
        ispending(P) = 1;
        actual(P) = i_above;
        dist(P) = d_above;
    end
    loop=loop+1;
    if loop>50
        break;
    end
end


% contract s2->b2 if needed:
if sum(G(b2,:)) == 1
    j = find(G(b2,:));
    G(s2,j) = 1;
    G_weight(s2,j) = G_weight(s2,b2) + G_weight(b2,j); 
    G(s2,b2) = 0; G_weight(s2,b2) = 0;
    G(b2,j) = 0; G_weight(b2,j) = 0;
end
% ignore isolated nodes (extra joining points):
needed = zeros(1,length(G)); % indicator of non-isolated node
for i=1:length(G)
    needed(i) = any([sum(G(i,:)) sum(G(:,i))]);
end
needed = logical(needed);
G = G(needed,needed);
G_weight = G_weight(needed,needed);

for i=1:length(G)
    for j=1:length(G)
        if G(i,j)>0 && G_weight(i,j)==0
            G_weight(i,j)=0.00001;
        end
    end
end


end

