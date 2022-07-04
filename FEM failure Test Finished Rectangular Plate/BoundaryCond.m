function [prescribedDof,activeDof] = BoundaryCond(nodeCoordinates,L,Wi,NBNd)

% Find Edges-Nodes
RIGHTNodes   = find(nodeCoordinates(:,1)==L)';
LEFTNodes    = find(nodeCoordinates(:,1)==0)';
UPNodes      = find(nodeCoordinates(:,2)==Wi)';
BELLOWtNodes = find(nodeCoordinates(:,2)==0)';
B = NBNd;

U  = unique(sort(LEFTNodes ));
V  = unique(sort(LEFTNodes+NBNd));

prescribedDof = [U V];
activeDof = setdiff([1:NBNd*2]' , prescribedDof);

end