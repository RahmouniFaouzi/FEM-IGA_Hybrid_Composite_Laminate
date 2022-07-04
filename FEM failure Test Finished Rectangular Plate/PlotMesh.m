function PlotMesh(coordinates,nodes,L,W)
%
% Plotting the Finite Element Mesh
% Initialization of the required matrices
nel  = size (nodes,1);
nnel = size (nodes,2);

X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;

% Extract X,Y coordinates for the (iel)-th element
for iel = 1:nel
    X(:,iel) = coordinates(nodes(iel,:),1) ;
    Y(:,iel) = coordinates(nodes(iel,:),2) ;
end
fh = figure;
set(fh,'name','Preprocessing for FEA','numbertitle','off','color','w') ;
patch(X,Y,'w')
title('Finite Element Mesh') ;
%axis off
axis([0. L*1.01 0. W*1.01])
if L==W
    axis equal ;
end