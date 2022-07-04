function [force] = formForceVectorKUV(L,nodeCoordinates,force,N,Ny)
% Linear application of forces my there is an erroe

RIGHTNodes   = find(nodeCoordinates(:,1)==L)';
LEFTNodes    = find(nodeCoordinates(:,1)==0)';


XG= [0.57735026918963D0, -0.57735026918963D0];
WGT=[1.00000000000000D0,  1.00000000000000D0];

if N(1)~= 0
    if N(2) == 1
        forceEdgeU = LEFTNodes  ;
    elseif N(2) == 2
        forceEdgeU = RIGHTNodes ;
    elseif N(2) == 3
        forceEdgeU = [LEFTNodes,RIGHTNodes] ;
    end
end

if N(1)== 0, return; end

if N(1)~= 0
    Y = forceEdgeU;
    j = 1;
    for ii = 1 : 2
        Ai = forceEdgeU(1:Ny+1);
        for i = 1 : length(Ai)-1
            if (ii == 2), i = i+ Ny +1; end
            frn(j,:) = [Y(i) Y(i+1)];
            j = j+1;
        end
    end
    
    for e=1:size(frn,1)
        if e > Ny, Aii = N(1); else Aii = -N(1);end
        sctr  = frn(e,:);
        sctrx = frn(e,:);
        for q=1:size(WGT,2)
            pt        = XG(q);
            wt        = WGT(q);
            [Ni,dNdxi]= L_basis(pt);
            J0        = dNdxi'*nodeCoordinates(sctr,:);
            detJ0     = norm(J0);
            force(sctrx) = force(sctrx)+Ni*Aii*detJ0*wt;
        end
    end
    
end
end
%%

function [N,dNdxi]= L_basis(coord)
%    1---------2
if size(coord,2) < 1
    disp('Error coordinate needed for the L2 element')
else
    xi=coord(1);
    N=([1-xi,1+xi]/2)';
    dNdxi=[-1;1]/2;
end
end

