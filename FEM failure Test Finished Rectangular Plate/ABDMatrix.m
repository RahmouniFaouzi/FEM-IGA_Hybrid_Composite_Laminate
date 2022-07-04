function[AMatrix,qbarra] =  ABDMatrix(Nplies, thetadt, h_ply, Prop,CLS)

Q = zeros(3,3,2);

for i = 1:2
    E1   = Prop(i,1); E2  = Prop(i,2);
    nu12 = Prop(i,3); G12 = Prop(i,4);
    nu21 = Prop(i,5);
    
    % Q matrix (material coordinates)
    denom = 1 - nu12 * nu21 ;
    Q11 = E1 / denom        ;
    Q12 = nu12 * E2 / denom ;
    Q22 = E2 / denom        ;
    Q66 = G12               ;
    
    Q(:,:,i) = [Q11  Q12  0;
                Q12  Q22  0;
                0    0  Q66];
end
% Qbar matrices (laminate coordinates) 
AMatrix = zeros(3,3);
qbarra  = zeros(3,3,Nplies);

for i = 1 : Nplies;
    m = cosd(thetadt(i)) ;
    n = sind(thetadt(i)) ;
    
    T = [ m^2  n^2   2*m*n; 
          n^2  m^2  -2*m*n; 
         -m*n  m*n  (m^2-n^2)];
     
    if CLS(i) == 1
        Qbar = inv(T) * Q(:,:,1) * (inv(T))' ; %#ok<*MINV>
    else 
        Qbar = inv(T) * Q(:,:,2) * (inv(T))' ;
    end
    
    qbarra(:,:,i) = Qbar;
    AMatrix = AMatrix + Qbar * h_ply;
end

end