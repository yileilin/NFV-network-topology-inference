function [ val ] = graph_edit_distance( G1, G2 )
% return the number of graph edits to change G1 into G2 under the optimal
% mapping of their vertices, subject to constraints that:
% known nodes (source/destination) are mapped exactly
% if a labeled node is mapped to a node of different label, an extra edit
% is incurred (label substitution). 
% Data structure: 
%   G.adjacency: adjacency matrix (sparse), 1: link, 0: no link 
%   G.source: 1*n array of source node indices; 
%   G.destination: 1*n array of destination node indices; 
%   G.vnf: 1*|F| cell array, G.vnf{i} are node indices with label f_i.
% Note: (i) G.vnf = [] if unlabeled (e.g., RNJ, CE); (ii) if G.vnf is not
% empty, length(G.vnf) must equal total #VNFs, nodes not in source/destination
% or any G.vnf{f} are considered dummies.
option = 1; %1: count each new node as one edit; 2: count each new node as two edits (one for labeling it)

if length(G1.adjacency) > length(G2.adjacency)
    tmp = G1; G1 = G2; G2 = tmp;
end% make sure G1 has fewer vertices than G2

if isempty(G1.vnf) || isempty(G2.vnf) % if either graph is unlabeled
    val = graph_edit_distance_unlabeled(G1, G2, option);
else % if both G1 and G2 are labeled
    val = graph_edit_distance_labeled1(G1, G2, option);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ val ] = graph_edit_distance_unlabeled(G1, G2, option)
% all nodes except for sources and destinations are mapped arbitrarily
% (brute-force search among possible mappings)
map = zeros(1,length(G1.adjacency)); % map(i): index of i-th node in G1 in G2
for i=1:length(G1.source)
    map(G1.source(i)) = G2.source(i);
    map(G1.destination(i)) = G2.destination(i); 
end% fix mapping for known nodes
% brute-force search of mapping of remaining nodes:
unmapped1 = find(map == 0); % unmapped nodes in G1
unmapped2 = setdiff([1:length(G2.adjacency)],[G2.source G2.destination]); % candidate nodes in G2
val = inf; 
allsets = allperms(unmapped2, length(unmapped1));
for i=1:length(allsets(:,1)) % for each subset of candidate nodes in each order
    map(unmapped1) = allsets(i,:);
    % augment G1.adjacency into a matrix of the same size as
    % G2.adjacency:
    G1aug = zeros(length(G2.adjacency));
    G1aug(map,map) = G1.adjacency; % G1aug is the new adjacency matrix of G1 after mapping nodes and adding new nodes
    edits = length(G2.adjacency) - length(G1.adjacency) + nnz(G1aug - G2.adjacency); % #edits for this mapping (including adding nodes, adding/removing links)
    if ~isempty(G1.vnf) || ~isempty(G2.vnf) % if one graph is labeled and the other is not
        edits = edits + length(unmapped1); % add #edits for node labels
    end
    if ~isempty(G2.vnf) && option == 2 % if the larger graph is labeled and the smaller graph is not
        edits = edits + length(G2.adjacency) - length(G1.adjacency); % add #edits for labeling newly added nodes
    end
    val = min(val, edits);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This implementation is slower but yields smaller edit distance, as it
% searches among all possible mappings of non-source/destination nodes. 
function [ val ] = graph_edit_distance_labeled1(G1, G2, option)
% allowing arbitrary mapping of non-source/destination nodes, but counting
% the edits of labels for mapped nodes
% (brute-force search among possible mappings)

% if (length(G1.adjacency)<length(G2.adjacency))
%     G=G1;
%     G1=G2;
%     G2=G;
% end
map = zeros(1,length(G1.adjacency)); % map(i): index of i-th node in G1 in G2
map(G1.source(1)) = G2.source(1);
for i=1:length(G1.destination)
    %map(G1.source(i)) = G2.source(i);
    map(G1.destination(i)) = G2.destination(i); 
end% fix mapping for known nodes
% brute-force search of mapping of remaining nodes:
unmapped1 = find(map == 0); % unmapped nodes in G1
unmapped2 = setdiff([1:length(G2.adjacency)],[G2.source G2.destination]); % candidate nodes in G2
label1 = find_label(G1); label2 = find_label(G2); % find label of each node in {1,...,F+1}
val = inf; 
MAX_ARRAY = 10^9; % maximum possible array (bytes) divided by 1 (by using 'uint8')
% if (length(unmapped2)<length(unmapped1))
%     unmapped3=unmapped2;
%     unmapped2=unmapped1;
%     unmapped1=unmapped3;
% end
if nchoosek(length(unmapped2),length(unmapped1))*factorial(length(unmapped1))*length(unmapped1) < MAX_ARRAY
    % brute-force search:
    allsets = allperms(unmapped2, length(unmapped1));
    for i=1:length(allsets(:,1)) % for each subset of candidate nodes in each order
        map(unmapped1) = allsets(i,:);
        % augment G1.adjacency into a matrix of the same size as
        % G2.adjacency:
        G1aug = zeros(length(G2.adjacency));
        G1aug(map,map) = G1.adjacency; % G1aug is the new adjacency matrix of G1 after mapping nodes and adding new nodes
        edits = length(G2.adjacency) - length(G1.adjacency) + nnz(G1aug - G2.adjacency); % #edits for this mapping (including adding nodes, adding/removing links)
        edits = edits + sum(label1(unmapped1) ~= label2(map(unmapped1))); % #edits for labels of mapped nodes
        if option == 2 % if the larger graph is labeled and the smaller graph is not
            edits = edits + length(G2.adjacency) - length(G1.adjacency); % add #edits for labeling newly added nodes
        end
        val = min(val, edits);
    end
else% too many mappings, random search:
    for i=1:10^6 %nchoosek(length(unmapped2),length(unmapped1))*factorial(length(unmapped1))        
        randset = randperm(length(unmapped2),length(unmapped1) ); randset = unmapped2(randset);
        map(unmapped1) = randset;
        % augment G1 into same size as G2:
        G1aug = zeros(length(G2.adjacency));
        G1aug(map,map) = G1.adjacency; % G1aug is the new adjacency matrix of G1 after mapping nodes and adding new nodes
        edits = length(G2.adjacency) - length(G1.adjacency) + nnz(G1aug - G2.adjacency); % #edits for this mapping (including adding nodes, adding/removing links)
        edits = edits + sum(label1(unmapped1) ~= label2(map(unmapped1))); % #edits for labels of mapped nodes
        if option == 2 % if the larger graph is labeled and the smaller graph is not
            edits = edits + length(G2.adjacency) - length(G1.adjacency); % add #edits for labeling newly added nodes
        end
        val = min(val, edits);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This implementation is faster but yields larger edit distance, as it only
% searches mappings that map nodes of the same label (before considering
% any mixed-label mapping). 
function [ val ] = graph_edit_distance_labeled(G1, G2, option)
isdebug = 0; % set to 1 for printing best mapping found 

% only map nodes of the same label
map = zeros(1,length(G1.adjacency)); % map(i): index of i-th node in G1 in G2
for i=1:length(G1.source)
    map(G1.source(i)) = G2.source(i);
    map(G1.destination(i)) = G2.destination(i); 
end% fix mapping for known nodes
unmapped1 = find(map == 0); % non-source/destination nodes
unmapped2 = setdiff([1:length(G2.adjacency)],[G2.source G2.destination]);
% brute-force search of mapping for remaining nodes:
val = inf;
F = length(G1.vnf); % #VNFs
allmaps = cell(1,F+1); % allmaps{f} contains all mapps of nodes with label f (including dummy)
for f=1:F % for each type of VNF
    if length(G1.vnf{f}) <= length(G2.vnf{f}) % map type-f nodes from G1 to G2
        allmaps{f} = allperms(G2.vnf{f}, length(G1.vnf{f})); % map(G1.vnf{f}) = allmaps{f}(i,:) 
    else % map type-f nodes from G2 to G1
        allmaps{f} = allperms(G1.vnf{f}, length(G2.vnf{f})); % map(allmaps{f}(i,:)) = G2.vnf{f}
    end    
end
G1.dummy = find_dummy(G1); G2.dummy = find_dummy(G2); % now find all mappings of dummies
if length(G1.dummy) <= length(G2.dummy) % map dummies from G1 to G2
    allmaps{F+1} = allperms(G2.dummy, length(G1.dummy)); % map(G1.dummy) = allmaps{F+1}(i,:)
else % map dummies from G2 to G1
    allmaps{F+1} = allperms(G1.dummy, length(G2.dummy)); % map(allmaps{F+1}(i,:)) = G2.dummy
end
n_maps = zeros(1,F+1); % n_maps(f): #mappings in category f
for f=1:F+1
    if ~isempty(allmaps{f})
        n_maps(f) = length(allmaps{f}(:,1));
    else
        n_maps(f) = 1; % only one choice: empty set
    end
end% #intra-category mappings
n_map = prod(n_maps); 
bestmap = []; 
for i=1:n_map % for each intra-category mapping
    I = intra_category_mapping(n_maps,i); % select allmaps{f}(I(f),:) for each f 
    is_unmapped1 = zeros(1,length(G1.adjacency)); is_unmapped1(unmapped1) = 1;
    is_unmapped2 = zeros(1,length(G2.adjacency)); is_unmapped2(unmapped2) = 1;
    for f=1:F
        if isempty(allmaps{f})
            continue;
        end
        if length(G1.vnf{f}) <= length(G2.vnf{f})
            map(G1.vnf{f}) = allmaps{f}(I(f),:);
            is_unmapped1(G1.vnf{f}) = 0;
            is_unmapped2(allmaps{f}(I(f),:)) = 0;
        else
            map(allmaps{f}(I(f),:)) = G2.vnf{f};
            is_unmapped1(allmaps{f}(I(f),:)) = 0;
            is_unmapped2(G2.vnf{f}) = 0;
        end
    end
    if length(G1.dummy) <= length(G2.dummy) % map dummies from G1 to G2
        map(G1.dummy) = allmaps{F+1}(I(F+1),:);
        is_unmapped1(G1.dummy) = 0;
        is_unmapped2(allmaps{F+1}(I(F+1),:)) = 0;
    else % map dummies from G2 to G1
        map(allmaps{F+1}(I(F+1),:)) = G2.dummy;
        is_unmapped1(allmaps{F+1}(I(F+1),:)) = 0;
        is_unmapped2(G2.dummy) = 0;
    end% done selecting intra-category mappings
    edits = sum(is_unmapped1) + (length(G2.adjacency) - length(G1.adjacency)); % #nodes edits, including changing labels and adding nodes
%     In option 2, each added node counts as 2 edits (add node, add its label).
    if option == 2
        edits = edits + (length(G2.adjacency) - length(G1.adjacency));
    end
    candidates = 1:length(G2.adjacency); candidates = candidates(is_unmapped2>0);
    allsets = allperms( candidates, sum(is_unmapped1) );
    if isempty(allsets)
        n_map_inter = 1;
    else
        n_map_inter = length(allsets(:,1));
    end
    for j=1:n_map_inter % for each inter-category mapping
        if ~isempty(allsets)
            map(is_unmapped1>0) = allsets(j,:);
        end
        % augment G1.adjacency into a matrix of the same size as
        % G2.adjacency:
        G1aug = zeros(length(G2.adjacency));
        G1aug(map,map) = G1.adjacency; % G1aug is the new adjacency matrix of G1 after mapping nodes and adding new nodes
        edits = edits + nnz(G1aug - G2.adjacency); % #link edits (including adding/removing links)
        if isdebug
            if edits < val
                bestmap = map;
            end
        end
        val = min(val, edits);
    end
end
if isdebug
    bestmap
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = allperms(v, n)
% all subsets of n items in v in all orders:
val = zeros(nchoosek(length(v),n)*factorial(n),n,'uint8');
if n>0
    allsets = nchoosek(v, n);
    for i=1:length(allsets(:,1))
        val((i-1)*factorial(n)+1: i*factorial(n),:) = uint8(perms(allsets(i,:)));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function dummy = find_dummy(G)
isdummy = ones(1,length(G.adjacency)); % isdummy(i)=1 iff node i is dummy
isdummy(G.source) = 0;
isdummy(G.destination) = 0;
for f=1:length(G.vnf)
    isdummy(G.vnf{f}) = 0;
end
dummy = find(isdummy);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = intra_category_mapping(n_maps,i)
% convert i-th combination to exact combination of mapping indices
I = zeros(1,length(n_maps));
divider = prod(n_maps);
for j=length(n_maps):-1:1
    divider = divider / n_maps(j); 
    I(j) = ceil(i/divider);
    i = i - divider*(I(j)-1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label = find_label(G)
% label(v) is label of node v: 0 for source/destination; {1,...,F+1} for
% other nodes
label = zeros(1,length(G.adjacency));
for f=1:length(G.vnf)
    label(G.vnf{f}) = f;
end
label(find_dummy(G)) = length(G.vnf)+1;
end