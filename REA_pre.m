function [tree, weight, R, pathlength, sharedlength, s2b2]=REA_pre(path_length,sharedpath_length)
weight=zeros(4,4);
weight(4,3)=sharedpath_length(1);
weight(3,1)=path_length(1)-sharedpath_length(1);
weight(3,2)=path_length(2)-sharedpath_length(1);
tree=zeros(4,4);
tree(4,3)=1;tree(3,1)=1;tree(3,2)=1;
pathlength=path_length(3:4);
sharedlength(1)=sharedpath_length(2);
sharedlength(2)=sharedpath_length(5);
R=[1,2];
s2b2=sharedpath_length(6);

%only consider two sources
% tree: adjacency matrix of a directed tree from s1 to all receivers R; tree(i,j)=1 means i->j edge
% exists
% weight: weight(i,j) is the weight of edge (i,j)
% R: node indices of receivers in tree
% pathlength: pathlength(i) is the length of s2-to-R(i) path
% sharedlength: sharedlength(i) is the shared length between path
% s1-to-R(i) and path s2-to-R(i)
% s2b2: length of link (s2,b2), i.e., minimum shared length between s2 and
% any two receivers in R
% global n_dest
% global n_source
% global n_path
% win=1:n_dest;
% path_length_one=path_length(win);
% temp=nchoosek(1:n_path,2);
% sharedpath_length_one=[];
% % for i=1:ceil(length(temp(:,1)))
% %     if ismember(temp(i,1),win)&&ismember(temp(i,2),win)
% %         sharedpath_length_one=[sharedpath_length_one,sharedpath_length(i)];
% %     end
% % end
% % [weight,~,R]=RNJ_plus(path_length_one,sharedpath_length_one);
% 
% tree=weight;
% % for i=1:length(tree)
% %     for j=1:length(tree)
% %         if tree(i,j)>0
% %             tree(i,j)=1;
% %         end
% %     end
% % end
% % pathlength=path_length((n_dest+1):(2*n_dest));
% % sharedlength=[];
% % for i=1:n_dest
% %     for j=1:length(temp(:,1))
% %         if (i==temp(j,1))&&(i+n_dest==temp(j,2))
% %             sharedlength=[sharedlength,sharedpath_length(j)];
% %         end
% %     end
% % end
% % win2=(n_dest+1):(2*n_dest);
% % shared2=[];
% % for i=1:length(temp(:,1))
% %     if ismember(temp(i,1),win2)&&ismember(temp(i,2),win2)
% %         shared2=[shared2,sharedpath_length(i)];
% %     end
% % end
% s2b2=min(shared2);

end