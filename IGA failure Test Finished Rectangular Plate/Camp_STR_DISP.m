function [stress, disp]  = Camp_STR_DISP(noElems , index, ...
    elRangeU , elRangeV ,element , noCtrPts, controlPts,noDofs, ...
    noPtsX , noPtsY , uKnot , vKnot ,weights , C , U, p , q)


elementV = build_Q4_elements( length(unique(uKnot)), length(unique(vKnot))) ;

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);
disp   = zeros(noElems,size(elementV,2),2);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    sctr   = element(e,:);
    pts    = controlPts(sctr,:);
    
    uspan = Find_span(noPtsX-1,p,xiE(1), uKnot);
    vspan = Find_span(noPtsY-1,q,etaE(1),vKnot);
    
    % loop over Gauss points
    
    gp = 1;
    for iv=1:2
        if (iv==2)
            xiE = sort(xiE,'descend');
        end
        for iu=1:2
            Xi  = xiE(iu);
            Eta = etaE(iv);
            
            [N, dRdxi, dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            jacob      = pts' * [dRdxi' dRdeta'];
            
            if (det(jacob) <= 1e-8)
                [N, ~, ~] = NURBS2DBasisDers([Xi-0.01; ...
                    Eta+0.001],p,q,uKnot,vKnot,weights');
            end
            
            disp(e,gp,:)    = N*[Ux(sctr) Uy(sctr)];
            
            gp = gp +1;
        end
    end
end
%
InD = size(C,3);
stress = zeros(noElems,size(elementV,2),3 , InD);

for M = 1:InD
    for e=1:noElems
        idu    = index(e,1);
        idv    = index(e,2);
        xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
        
        sctr   = element(e,:);         %  element scatter vector
        sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
        nn     = length(sctr);
        
        B      = zeros(3,2*nn);
        pts    = controlPts(sctr,:);
        
        uspan = Find_span(noPtsX-1,p,xiE(1), uKnot);
        vspan = Find_span(noPtsY-1,q,etaE(1),vKnot);
        
        % loop over Gauss points
        
        gp = 1;
        for iv=1:2
            if (iv==2)
                xiE = sort(xiE,'descend');
            end
            for iu=1:2
                Xi  = xiE(iu);
                Eta = etaE(iv);
                [~, dRdxi, dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                    p,q,uKnot,vKnot,weights',[uspan;vspan]);
                
                % compute the jacobian of physical and parameter domain mapping
                % then the derivative w.r.t spatial physical coordinates
                
                jacob      = pts' * [dRdxi' dRdeta'];
                
                if (det(jacob) <= 1e-8)
                    [~, dRdxi, dRdeta] = NURBS2DBasisDers([Xi-0.01; ...
                        Eta+0.001],p,q,uKnot,vKnot,weights');
                    jacob      = pts' * [dRdxi' dRdeta'];
                end
                
                % Jacobian inverse and spatial derivatives
                
                dRdx       = [dRdxi' dRdeta']/jacob;
                
                % B matrix
                
                B(1,1:nn)       = dRdx(:,1)';
                B(2,nn+1:2*nn)  = dRdx(:,2)';
                B(3,1:nn)       = dRdx(:,2)';
                B(3,nn+1:2*nn)  = dRdx(:,1)';
                
                strain          = B*U(sctrB);
                
                stress(e,gp,:,M)  = C(:,:,M)*strain;
                
                gp = gp +1;
            end
        end
    end
end
%
end


%

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
%
function knotSpanIndex = Find_span(n,p,u,U)
if (u == U(n+2))
    knotSpanIndex= n;
    return
end
low = p;
high = n+1;
mid = floor((low + high)/2);
while (u <U(mid+1) || u >= U(mid+2) )
    if( u < U(mid+1))
        high = mid;
    else
        low = mid;
    end
    mid = floor((low+high)/2);
end
knotSpanIndex = mid;
end