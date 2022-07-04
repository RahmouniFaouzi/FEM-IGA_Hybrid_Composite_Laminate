function Bm = call_Bm( xi,eta,elemType, nodeCoordinates, ndof)

[~ ,naturalDerivatives]  =  shapeFunctionsQ(xi,eta,elemType);
[ ~ , ~ , XYderivatives] =  Jacobian(nodeCoordinates,naturalDerivatives);

Bm = zeros (3,ndof); 
Bm(1,1:4) = XYderivatives(:,1)'; 
Bm(2,5:8) = XYderivatives(:,2)'; 
Bm(3,1:4) = XYderivatives(:,2)'; 
Bm(3,5:8) = XYderivatives(:,1)'; 

end