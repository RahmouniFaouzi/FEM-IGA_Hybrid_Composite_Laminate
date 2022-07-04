function [node ,elementV] = buildVisualizationMesh ( controlPts, weights,...
                         uKnot, vKnot, noPtsX, noPtsY , p, q )

% get rid of zero measure knot spans
uKnotVec = unique(uKnot);
vKnotVec = unique(vKnot);

% number of distinct knot values
noKnotsU = length(uKnotVec);
noKnotsV = length(vKnotVec);

elementV = build_Q4_elements(noKnotsU,noKnotsV);

% using homogeneous coordinates
projcoord = nurb2proj(noPtsX*noPtsY, controlPts, weights');

dim=size(projcoord,2);

node  = zeros(noKnotsU*noKnotsV,2);
count = 1;

for vk=1:noKnotsV
    eta = vKnotVec(vk);
    for uk=1:noKnotsU
        xi = uKnotVec(uk);
        tem = SurfacePoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
                           projcoord,dim,xi,eta);
        node(count,1) = tem(1)/tem(3);
        node(count,2) = tem(2)/tem(3);
        count = count + 1;
    end
end

end

function elementV = build_Q4_elements(noKnotsU,noKnotsV)
nnx   = noKnotsU;
nny   = noKnotsV;
inc_u = 1;
inc_v = nnx;
node_pattern = [1 2 nnx+2 nnx+1];
elementV     = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v);

end
%
function element = make_elem(node_pattern,num_u,num_v,inc_u,inc_v)
% creates a connectivity list

if ( nargin < 5 )
   disp(['Not enough parameters specified for make_elem function'])
end

inc=[zeros(1,size(node_pattern,2))];
e=1;
element=zeros(num_u*num_v,size(node_pattern,2));

for row=1:num_v
   for col=1:num_u
      element(e,:)=node_pattern+inc;
      inc=inc+inc_u;
      e=e+1;
   end
   inc=row*inc_v;
end
end